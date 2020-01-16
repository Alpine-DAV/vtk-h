#include <vtkh/filters/LagrangianInterpolation.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include <sstream>
#include <vector>
#include <cmath>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/Integrators.h>
#include <vtkm/worklet/particleadvection/Particles.h>


#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif

namespace vtkh 
{

LagrangianInterpolation::LagrangianInterpolation()
{
}

LagrangianInterpolation::~LagrangianInterpolation()
{
}

void 
LagrangianInterpolation::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

//void 
//LagrangianInterpolation::SetInputPath(const std::string &input_path)
//{
//  m_input_path = input_path;
//}

void 
LagrangianInterpolation::SetWriteFrequency(const int &write_frequency)
{
  m_write_frequency = write_frequency;
}

void LagrangianInterpolation::PreExecute() 
{
  Filter::PreExecute();
}

void LagrangianInterpolation::PostExecute()
{
  Filter::PostExecute();
}

void LagrangianInterpolation::DoExecute()
{
/*
* Algorithm
* - Get MPI rank, number of processes
* - Based on input parameters: input path, write frequency -> load corresponding VTK file
* If number of processes greater than 1:
* * - Extract BBOX of node specific VTK file
* * - Create an array to store all BBOX of all VTK files
* * - Share BBOX information with all processes
* * - Calculate a set of query points (26) to identify neighbors 
* * - Evaluate query points against BBOX information to create an adjacent list
* * - Load all adjacent VTK files and own VTK file to create a basis flow list
* * - Nested loop - loop over all points to find invalid points
* * - For each invalid point find all points within a distance X to use for shepard's interpolation.
* * - Change value of each "displacement" field value where "valid" is 0 -> Use Shepard's interpolation. 
* Load seed array.
* Identify local seeds.
* Perform VTK-m particle advection
* Exchange particles as necessary. Loop continues. 
*/

#ifdef VTKH_PARALLEL

  vtkm::Id rank = vtkh::GetMPIRank();
  vtkm::Id num_ranks = vtkh::GetMPISize();
  bool allReceived[num_ranks] = {false};
  allReceived[rank] = true;

  float radius = 0.2;  // TODO User parameter radius
  int num_seeds = 1000; // TODO User parameter total number of seeds

  std::vector<double> px, py, pz; // Start location of a basis flow
  std::vector<double> dx, dy, dz; // Displacement of corresponding basis flow

  std::string seed_input_path = "/research/Sudhanshu/Cloverleaf3D/Seeds.txt"; // TODO User parameter
  std::ifstream seed_stream(seed_input_path);
  
  float x1, y1, z1;
  int seed_counter = 0;
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>> SeedParticleArray, SeedParticleOriginal;
  SeedParticleArray.Allocate(num_seeds);
  SeedParticleOriginal.Allocate(num_seeds);
  auto seed_portal = SeedParticleArray.GetPortalControl();
  std::vector<int> seed_validity;
  std::vector<int> seed_rank;
  std::vector<vtkm::Id> seed_id;

  while(seed_stream >> x1)
  {
    seed_stream >> y1;
    seed_stream >> z1;
    seed_portal.Set(seed_counter, vtkm::Vec<vtkm::Float64, 3>(x1,y1,z1));
    seed_id.push_back(seed_counter);
    seed_counter++;
  }

  vtkm::cont::ArrayCopy(SeedParticleArray, SeedParticleOriginal);
  auto seed_original_portal = SeedParticleOriginal.GetPortalControl();

  std::string input_path = "/research/Sudhanshu/Cloverleaf3D/"; // TODO User parameter
  int interval = 10; // TODO User parameter write_frequency
  int start_cycle = 10; // TODO User parameter start
  int end_cycle = 200; // TODO User parameter end
  int cycle = interval;
  
  std::vector<int> neighbor_ranks;
    
  double *bbox_list = (double*)malloc(sizeof(double)*6*num_ranks);

  for(int cycle = start_cycle; cycle <= end_cycle; cycle += interval)
  {
    if(rank == 0)
      std::cout << "[" << rank << "] Starting cycle: " << cycle << std::endl;
    std::stringstream filename;
    filename << input_path << "Lagrangian_flowmap_" << rank << "_" << cycle << ".vtk";
    
    vtkm::cont::DataSet input_flowmap;
    vtkm::io::reader::VTKDataSetReader reader(filename.str().c_str());
    input_flowmap = reader.ReadDataSet();

    int num_pts = input_flowmap.GetNumberOfPoints();
    auto coords_portal = input_flowmap.GetCoordinateSystem().GetData().GetPortalControl();      

    auto valid_VAH = input_flowmap.GetField("valid").GetData();
    vtkm::cont::ArrayHandle<vtkm::Int32> valid_arrayhandle;
    valid_VAH.CopyTo(valid_arrayhandle);

    auto disp_VAH = input_flowmap.GetField("displacement").GetData();
    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64, 3>> disp_arrayhandle;
    disp_VAH.CopyTo(disp_arrayhandle);

    auto valid_portal = valid_arrayhandle.GetPortalControl();
    auto disp_portal = disp_arrayhandle.GetPortalControl();      

    vtkm::Bounds bounds = input_flowmap.GetCoordinateSystem().GetBounds();
    double BB[6];
    BB[0] = bounds.X.Min;  
    BB[1] = bounds.X.Max;  
    BB[2] = bounds.Y.Min;  
    BB[3] = bounds.Y.Max;  
    BB[4] = bounds.Z.Min;  
    BB[5] = bounds.Z.Max;  
      
    

    // Each rank needs to identify if a seed belongs to its node
    if(cycle == start_cycle)
    { // Initialize seed validity on the first interval
      for(int i = 0; i < num_seeds; i++)
      {
        auto pt = seed_portal.Get(i);
        if(BoundsCheck(pt[0], pt[1], pt[2], BB))
        {
          seed_validity.push_back(1);
        }
        else
        {
          seed_validity.push_back(0);
        }
      }
    }
      
    if(num_ranks > 1)
    {
      if(cycle == start_cycle)
      { // Initialize seed validity on the first interval
        bbox_list[rank*6 + 0] = bounds.X.Min;
        bbox_list[rank*6 + 1] = bounds.X.Max;
        bbox_list[rank*6 + 2] = bounds.Y.Min;
        bbox_list[rank*6 + 3] = bounds.Y.Max;
        bbox_list[rank*6 + 4] = bounds.Z.Min;
        bbox_list[rank*6 + 5] = bounds.Z.Max;
        
        for(int i = 0; i < num_ranks; i++)
        {
          if(i != rank)
          {
            int ierr = MPI_Bsend(BB, 6, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
          }
        }

        while(!AllMessagesReceived(allReceived, num_ranks))
        {
          MPI_Status probe_status, recv_status;
          int ierr = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &probe_status);
          int count;
          MPI_Get_count(&probe_status, MPI_DOUBLE, &count);
          double *recvbuff;
          recvbuff = (double*)malloc(sizeof(double)*count);
          MPI_Recv(recvbuff, count, MPI_DOUBLE, probe_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_status);
          if(count == 6)
          {
            bbox_list[6*probe_status.MPI_SOURCE + 0] = recvbuff[0];
            bbox_list[6*probe_status.MPI_SOURCE + 1] = recvbuff[1];
            bbox_list[6*probe_status.MPI_SOURCE + 2] = recvbuff[2];
            bbox_list[6*probe_status.MPI_SOURCE + 3] = recvbuff[3];
            bbox_list[6*probe_status.MPI_SOURCE + 4] = recvbuff[4];
            bbox_list[6*probe_status.MPI_SOURCE + 5] = recvbuff[5];

            allReceived[recv_status.MPI_SOURCE] = true;
          }
          else
          {
            std::cout << "[" << rank << "] Corrupt message received from " << probe_status.MPI_SOURCE << std::endl;
          }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // Calculate query points and for list of adjacent nodes  
        neighbor_ranks = GetNeighborRankList(num_ranks, rank, bbox_list);  
      } 

     
      /* Scan input_flowmap */
      for(int i = 0; i < num_pts; i++)
      {
        auto pt = coords_portal.Get(i);
        auto disp = disp_portal.Get(i);
        if(valid_portal.Get(i)) // If valid == 1
        {
          px.push_back(pt[0]);
          py.push_back(pt[1]);
          pz.push_back(pt[2]);

          dx.push_back(disp[0]);
          dy.push_back(disp[1]);
          dz.push_back(disp[2]);
        } 
      } // Loop over points of the node-specific input flowmap.  

      /* Loop over neighbor_ranks */
      // VTK objects to read flowmaps of neighboring ranks
      vtkm::cont::DataSet neighbor_flowmap;
      for(int n = 0; n < neighbor_ranks.size(); n++)
      {
        std::stringstream neighbor_filename;
        neighbor_filename << input_path << "Lagrangian_flowmap_" << neighbor_ranks[n] << "_" << interval << ".vtk";
    
        vtkm::io::reader::VTKDataSetReader neighbor_reader(neighbor_filename.str().c_str());
        neighbor_flowmap = neighbor_reader.ReadDataSet();

        int neighbor_pts = neighbor_flowmap.GetNumberOfPoints();
        auto n_coords_portal = neighbor_flowmap.GetCoordinateSystem().GetData().GetPortalControl();      
        auto n_valid_VAH = neighbor_flowmap.GetField("valid").GetData();
        vtkm::cont::ArrayHandle<vtkm::Int32> n_valid_arrayhandle;
        n_valid_VAH.CopyTo(n_valid_arrayhandle);
 
        auto n_disp_VAH = neighbor_flowmap.GetField("displacement").GetData();
        vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64, 3>> n_disp_arrayhandle;
        n_disp_VAH.CopyTo(n_disp_arrayhandle);
 
        auto n_valid_portal = n_valid_arrayhandle.GetPortalControl();
        auto n_disp_portal = n_disp_arrayhandle.GetPortalControl();      
 
        for(int i = 0; i < neighbor_pts; i++)
        {
          auto pt = n_coords_portal.Get(i);
          auto disp = n_disp_portal.Get(i);
          if(n_valid_portal.Get(i)) // If valid == 1
          {
            px.push_back(pt[0]);
            py.push_back(pt[1]);
            pz.push_back(pt[2]);

            dx.push_back(disp[0]);
            dy.push_back(disp[1]);
            dz.push_back(disp[2]);
          } 
        } // Loop over points
      } // Loop over neighboring ranks
    
    // Reconstruction is only needed when there are multiple nodes. Validity is consequential. 

    // We could probably use worklets to perform this operation in parallel. 
      int num_invalid = 0;
      for(int i = 0; i < num_pts; i++)
      {
        if(valid_portal.Get(i) == 0)
        {
          num_invalid++;
          auto query_pt = coords_portal.Get(i);
          std::vector<float> weights;
          float sum_weight, wx, wy, wz;
          wx = 0.0;
          wy = 0.0;
          wz = 0.0;
          int zeroDistIndex = -1;

          for(int j = 0; j < px.size(); j++)
          {
            float dist = sqrt(pow(query_pt[0] - px[j],2) + pow(query_pt[1] - py[j],2) + pow(query_pt[2] - pz[j],2));
            if(dist == 0.0)
            {
              zeroDistIndex = j;
              break;
            }
            float w = pow((std::max(0.0f, (radius - dist))/(radius*dist)),2);
            weights.push_back(w);
            sum_weight += w;
          }
            
          auto d = disp_portal.Get(i);
        
          if(zeroDistIndex >= 0)
          {
            std::cout << "[" << rank << "] ZDI Displacement at " << i << " was : " << d[0] << " " << d[1] << " " << d[2] << 
            " at " << query_pt[0] << " " << query_pt[1] << " " << query_pt[2] << " and is now : " <<
            dx[zeroDistIndex] << " " << dy[zeroDistIndex] << " " << dz[zeroDistIndex] << std::endl;

            disp_portal.Set(i, vtkm::Vec<vtkm::Float64, 3>(dx[zeroDistIndex], dy[zeroDistIndex], dz[zeroDistIndex])); 
          }
          else
          {
            for(int j = 0; j < px.size(); j++)
            {
              wx += (weights[j]*dx[j])/sum_weight;
              wy += (weights[j]*dy[j])/sum_weight;
              wz += (weights[j]*dz[j])/sum_weight;
            }
            std::cout << "[" << rank << "] NonZDI Displacement at " << i << " was : " << d[0] << " " << d[1] << " " << d[2] << 
            " at " << query_pt[0] << " " << query_pt[1] << " " << query_pt[2] << " and is now : " <<
            wx << " " << wy << " " << wz << std::endl;
            disp_portal.Set(i, vtkm::Vec<vtkm::Float64, 3>(wx, wy, wz)); 
          } 
        }
      }
    } // endif there are more than 1 node. 
    else{} // If there is only a single node -> then use the displacement field as is without any reconstruction.
      MPI_Barrier(MPI_COMM_WORLD);
      if(rank == 0)
        std::cout << "[" << rank << "] Reconstruction completed for all nodes" << std::endl;

    /* Particle Advection Start */

      using FieldHandle = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64, 3>>;
      const vtkm::cont::DynamicCellSet& cells = input_flowmap.GetCellSet();
      const vtkm::cont::CoordinateSystem& coords = input_flowmap.GetCoordinateSystem();

      using GridEvalType = vtkm::worklet::particleadvection::GridEvaluator<FieldHandle>;

      vtkm::worklet::ParticleAdvection particleadvection;
      vtkm::worklet::ParticleAdvectionResult res, ex_res;

      using EulerType = vtkm::worklet::particleadvection::EulerIntegrator<GridEvalType>;
      GridEvalType eval(coords, cells, disp_arrayhandle);
      EulerType euler(eval, static_cast<vtkm::Float32>(1)); // step size set to 1. Lagrangian-based advection.
      res = particleadvection.Run(euler, SeedParticleArray, 1);
 
    if(rank == 0)
      std::cout << "[" << rank << "] Executed particle advection worklet" << std::endl;
      MPI_Barrier(MPI_COMM_WORLD);

      auto particle_positions = res.positions;
      auto particle_stepstaken = res.stepsTaken;

      auto position_portal = particle_positions.GetPortalControl();
      auto steps_portal = particle_stepstaken.GetPortalControl();

    /* Particle Advection End */
    /* Write Output Start */

      int connectivity_index = 0;
      std::vector<vtkm::Id> connectivity;
      std::vector<vtkm::Vec<vtkm::Float64, 3>> pointCoordinates;
      std::vector<vtkm::UInt8> shapes;
      std::vector<vtkm::IdComponent> numIndices;
      std::vector<vtkm::Id> pathlineId;

      for(int i = 0; i < num_seeds; i++)
      {
        if(seed_validity[i])
        {
          auto pt1 = seed_original_portal.Get(i);
          auto pt2 = position_portal.Get(i);

          connectivity.push_back(connectivity_index);
          connectivity.push_back(connectivity_index + 1);
          connectivity_index += 2;
          pointCoordinates.push_back(
            vtkm::Vec<vtkm::Float64, 3>(pt1[0], pt1[1], pt1[2]));
          pointCoordinates.push_back(
            vtkm::Vec<vtkm::Float64, 3>(pt2[0], pt2[1], pt2[2]));
          shapes.push_back(vtkm::CELL_SHAPE_LINE);
          numIndices.push_back(2);
          pathlineId.push_back(i);
        }
      }
      
      vtkm::cont::DataSetBuilderExplicit DSB_Explicit;
      vtkm::cont::DataSet pathlines = DSB_Explicit.Create(pointCoordinates, shapes, numIndices, connectivity); 
      vtkm::cont::DataSetFieldAdd dsfa;
      dsfa.AddCellField(pathlines, "ID", pathlineId);

      std::stringstream outputfile;
      outputfile << "/research/Sudhanshu/Cloverleaf3D/Pathlines/Pathlines_" << rank << "_" << cycle << ".vtk";
      
      vtkm::io::writer::VTKDataSetWriter writer(outputfile.str().c_str());
      writer.WriteDataSet(pathlines);
    if(rank == 0)
      std::cout << "[" << rank << "] Wrote data to disk " << std::endl;
    /* Write Output End */
      
    /* Exchange information - Update Seed Information Start */

    std::vector<int> outgoing_id;
    std::vector<int> outgoing_dest;

    for(int i = 0; i < num_seeds; i++)
    {
      if(seed_validity[i])
      {
        auto pt = seed_portal.Get(i);
        if(!BoundsCheck(pt[0], pt[1], pt[2], BB))
        {
          seed_validity[i] = 0;
          int dest = -1;
          for(int r = 0; r < num_ranks; r++)
          {
            if(r != rank)
            {
              if(BoundsCheck(pt[0], pt[1], pt[2],
           bbox_list[r*6 + 0], bbox_list[r*6 + 1], bbox_list[r*6 + 2], bbox_list[r*6 + 3], bbox_list[r*6 + 4], bbox_list[r*6 + 5]))
              {
                dest = r;
                break;
              }
            }
          }
          if(dest != -1)
          {
            outgoing_id.push_back(i);
            outgoing_dest.push_back(dest);    
          }
        }
      }
    }

    // Send Receive particles
    if(num_ranks > 1)
    {
      MPI_Barrier(MPI_COMM_WORLD);

      int bufsize = 10*(num_seeds*5 + (MPI_BSEND_OVERHEAD * num_ranks)); // message packet size 4 + message packet size 1 for empty sends.
      double *buf = (double*)malloc(sizeof(double)*bufsize);
      MPI_Buffer_attach(buf, bufsize);  
      // Sending messages to all ranks.
      // Message format: ID X Y Z ID X Y Z ..
      // If not particles to send: -1 
      for(int r = 0; r < num_ranks; r++)
      {
        if(r != rank)
        {
          int particles_to_send = 0;
          std::vector<double> message;
          for(int j = 0; j < outgoing_dest.size(); j++)
          {
            if(outgoing_dest[j] == r)
            {
              auto pt = seed_portal.Get(outgoing_id[j]);
              particles_to_send++;
              message.push_back(outgoing_id[j]);
              message.push_back(pt[0]);
              message.push_back(pt[1]);
              message.push_back(pt[2]);
            }
          }
          
          int buffsize = particles_to_send*4;
          double *sendbuff;
          
          
          if(buffsize == 0)
            sendbuff = (double*)malloc(sizeof(double));
          else
            sendbuff = (double*)malloc(sizeof(double)*buffsize); 

          if(buffsize  == 0)
          {
            buffsize = 1;
            sendbuff[0] = -1.0;
          }
          else
          {
            for(int k = 0; k < message.size(); k++)
            {
              sendbuff[k] = message[k];
            }
          }
          int ierr = MPI_Bsend(sendbuff, buffsize, MPI_DOUBLE, r, 13, MPI_COMM_WORLD);
          free(sendbuff);
        } 
      }   // Loop over all ranks.
      MPI_Barrier(MPI_COMM_WORLD);

      // Receive Messages
      bool allReceived2[num_ranks] = {false};
      allReceived2[rank] = true;

      while(!AllMessagesReceived(allReceived2, num_ranks))
      {
        MPI_Status probe_status, recv_status;
        int ierr = MPI_Probe(MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, &probe_status);
        int count;
        MPI_Get_count(&probe_status, MPI_DOUBLE, &count);
        double *recvbuff;
        recvbuff = (double*)malloc(sizeof(double)*count);
        MPI_Recv(recvbuff, count, MPI_DOUBLE, probe_status.MPI_SOURCE, 13, MPI_COMM_WORLD, &recv_status);
        if(count == 1)
        {
          allReceived2[recv_status.MPI_SOURCE] = true;
        }
        else if(count % 4 == 0)
        {
          int num_particles = count/4;
          
          for(int i = 0; i < num_particles; i++)
          {
            int index = recvbuff[i*4+0];
            seed_validity[index] = 1;
            seed_portal.Set(index, vtkm::Vec<vtkm::Float64, 3>(recvbuff[i*4+1], recvbuff[i*4+2], recvbuff[i*4+3]));
          }
          allReceived2[recv_status.MPI_SOURCE] = true;
        }
        else
        {
          std::cout << "[" << rank << "] Received message of invalid length from : " << recv_status.MPI_SOURCE << std::endl;
          allReceived2[recv_status.MPI_SOURCE] = true;
        } 
        free(recvbuff);
      }
      MPI_Buffer_detach(&buf, &bufsize);
      free(buf);
      MPI_Barrier(MPI_COMM_WORLD);
    } // endif rank > 1
    
    vtkm::cont::ArrayCopy(SeedParticleArray, SeedParticleOriginal);
  /* Exchange information - Update Seed Information End */
  } // End loop over interval


  this->m_output = new DataSet();

  const int num_domains = this->m_input->GetNumberOfDomains();

  vtkh::DataSet ds;
  
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    // insert interesting stuff
    m_output->AddDomain(dom, domain_id);
  }
 /* 
  vtkh::LagrangianReconstruction l_reconstruct;
  l_reconstruct.SetField("braid");
  l_reconstruct.SetInput(m_output);
  l_reconstruct.Update();
  vtkh::DataSet *reconstructed_ds = l_reconstruct.GetOutput();
*/
#endif
}

inline bool LagrangianInterpolation::BoundsCheck(float x, float y, float z, double *BB)
{
  if(x >= BB[0] && x <= BB[1] && y >= BB[2] && y <= BB[3] && z >= BB[4] && z <= BB[5])
  {
    return true;
  }
  return false;
}

inline bool LagrangianInterpolation::BoundsCheck(float x, float y, float z, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
  if(x >= xmin && x <= xmax && y >= ymin && y <= ymax && z >= zmin && z <= zmax)
  {
    return true;
  }
  return false;
}


inline bool LagrangianInterpolation::AllMessagesReceived(bool *a, int num_ranks)
{
for(int i = 0; i < num_ranks; i++)
{
  if(a[i] == false)
    return false;
}
return true;
}

std::vector<int> LagrangianInterpolation::GetNeighborRankList(vtkm::Id num_ranks, vtkm::Id rank, double *bbox_list)
{
    /*  
    * For each node. Calculate 26 query points. Each of these query points is outside the node. 
    * Check if that query point exists in any of the other nodes. 
    * If none, then do nothing. If yes, then add the rank to the vector.
    */
    
  double xlen = bbox_list[rank*6+1] - bbox_list[rank*6+0];
  double ylen = bbox_list[rank*6+3] - bbox_list[rank*6+2];
  double zlen = bbox_list[rank*6+5] - bbox_list[rank*6+4];
  
  double mid[3];
  mid[0] = bbox_list[rank*6+0] + (xlen/2.0);
  mid[1] = bbox_list[rank*6+2] + (ylen/2.0);
  mid[2] = bbox_list[rank*6+4] + (zlen/2.0);
    
  double query_points[26][3] = {{mid[0] + xlen, mid[1], mid[2]},
                                {mid[0] - xlen, mid[1], mid[2]},
                                {mid[0] + xlen, mid[1] + ylen, mid[2]},
                                {mid[0] + xlen, mid[1] - ylen, mid[2]},
                                {mid[0] - xlen, mid[1] + ylen, mid[2]},
                                {mid[0] - xlen, mid[1] - ylen, mid[2]},
                                {mid[0], mid[1] + ylen, mid[2]},
                                {mid[0], mid[1] - ylen, mid[2]}, 
                                {mid[0] + xlen, mid[1], mid[2] + zlen},
                                {mid[0] - xlen, mid[1], mid[2] + zlen}, 
                                {mid[0] + xlen, mid[1] + ylen, mid[2] + zlen},
                                {mid[0] + xlen, mid[1] - ylen, mid[2] + zlen},
                                {mid[0] - xlen, mid[1] + ylen, mid[2] + zlen},
                                {mid[0] - xlen, mid[1] - ylen, mid[2] + zlen},
                                {mid[0], mid[1] + ylen, mid[2] + zlen},
                                {mid[0], mid[1] - ylen, mid[2] + zlen},
                                {mid[0], mid[1], mid[2] + zlen}, 
                                {mid[0] + xlen, mid[1], mid[2] - zlen},
                                {mid[0] - xlen, mid[1], mid[2] - zlen}, 
                                {mid[0] + xlen, mid[1] + ylen, mid[2] - zlen},
                                {mid[0] + xlen, mid[1] - ylen, mid[2] - zlen},
                                {mid[0] - xlen, mid[1] + ylen, mid[2] - zlen},
                                {mid[0] - xlen, mid[1] - ylen, mid[2] - zlen},
                                {mid[0], mid[1] + ylen, mid[2] - zlen},
                                {mid[0], mid[1] - ylen, mid[2] - zlen},
                                {mid[0], mid[1], mid[2] - zlen}};

  std::vector<int> neighbor_ranks;
  
  for(int q = 0; q < 26; q++)
  {   
    for(int p = 0; p < num_ranks; p++)
    {   
      if(p != rank)
      {   
        if(BoundsCheck(query_points[q][0], query_points[q][1], query_points[q][2],
        bbox_list[p*6+0], bbox_list[p*6+1], bbox_list[p*6+2], bbox_list[p*6+3], bbox_list[p*6+4], bbox_list[p*6+5]))
        {   
          neighbor_ranks.push_back(p);
          break;
        }
      }
    }
  }
  
  return neighbor_ranks;
}

std::string
LagrangianInterpolation::GetName() const
{
  return "vtkh::LagrangianInterpolation";
}

} //  namespace vtkh
