#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <vtkh/vtkh.hpp>
#include <vtkh/Error.hpp>
#include <vtkh/filters/LagrangianInterpolation.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkm/Types.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandleVirtualCoordinates.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/worklet/particleadvection/Integrators.h>
#include <vtkm/worklet/particleadvection/Particles.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandleCast.h>
#include <vtkm/cont/Invoker.h>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif

namespace vtkh
{

namespace worklets
{
  class ShepardInterpolation : public vtkm::worklet::WorkletMapField
  {
    public:
      using ControlSignature = void(FieldIn, FieldIn, FieldInOut);
      using ExecutionSignature = void(_1, _2, _3);
      
      // Some variables I will want to pass in the constructor. 
      
      VTKM_CONT ShepardInterpolation(double r, int n, std::vector<double> pX, std::vector<double> pY, std::vector<double> pZ, 
        std::vector<double> dX, std::vector<double> dY, std::vector<double> dZ)
      {
        radius = r;
        num_basis = n;
        px = pX;
        py = pY;
        pz = pZ;
        dx = dX;
        dy = dY;
        dz = dZ;
      }
  
      template <typename A, typename T>
      VTKM_EXEC void operator()(const A &validity, const T &coordinates, vtkm::Vec<vtkm::Float64, 3> &displacement) const
      {
        if(validity == 0)
        {
          double query_pt[3];
          query_pt[0] = coordinates[0];
          query_pt[1] = coordinates[1];
          query_pt[2] = coordinates[2];

          std::vector<float> weights;
          float sum_weight, wx, wy, wz;
          wx = 0.0;
          wy = 0.0;
          wz = 0.0;
          int zeroDistIndex = -1;
  
          for(int j = 0; j < num_basis; j++)
          {
            float dist = sqrt(pow(query_pt[0] - px[j],2) + pow(query_pt[1] - py[j],2) + pow(query_pt[2] - pz[j],2));
            if(dist == 0.0)
            {
              zeroDistIndex = j;
              break;
            }
            float w = pow((std::max(0.0d, (radius - dist))/(radius*dist)),2);
            weights.push_back(w);
            sum_weight += w;
          }
  
          if(zeroDistIndex >= 0)
          {
             displacement[0] = dx[zeroDistIndex]; 
             displacement[1] = dy[zeroDistIndex]; 
             displacement[2] = dz[zeroDistIndex]; 
          }
          else
          {
            for(int j = 0; j < px.size(); j++)
            {
              wx += (weights[j]*dx[j])/sum_weight;
              wy += (weights[j]*dy[j])/sum_weight;
              wz += (weights[j]*dz[j])/sum_weight;
            }
            displacement[0] = wx; 
            displacement[1] = wy;
            displacement[2] = wz;
          }
        }
        else
        {
          displacement[0] = displacement[0];
          displacement[1] = displacement[1];
          displacement[2] = displacement[2];
        }
      }

      private:
      double radius;
      int num_basis;
      std::vector<double> px, py, pz, dx, dy, dz;
  };

  class SeedValidityCheck : public vtkm::worklet::WorkletMapField
  {
    public:
      using ControlSignature = void(FieldIn, FieldOut);
      using ExecutionSignature = void(_1, _2);
      
      VTKM_CONT SeedValidityCheck(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
      {
        BBox[0] = xmin;
        BBox[1] = xmax;
        BBox[2] = ymin;
        BBox[3] = ymax;
        BBox[4] = zmin;
        BBox[5] = zmax;
      }
  
      template <typename A, typename T>
      VTKM_EXEC void operator()(const A &particle, T &validity) const
      {
        auto p = particle.Pos;
        
        if(p[0] >= BBox[0] && p[0] <= BBox[1] && p[1] >= BBox[2] && p[1] <= BBox[3] && p[2] >= BBox[4] && p[2] <= BBox[5])
        {
          validity = 1;
        }
        else
        {
          validity = 0;
        } 
      }

    private:
      double BBox[6]; 
  };
}

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

void LagrangianInterpolation::SetSeedPath(const std::string &seed_path)
{
  m_seed_path = seed_path;
}

void LagrangianInterpolation::SetOutputPath(const std::string &output_path)
{
  m_output_path = output_path;
}

void LagrangianInterpolation::SetBasisPath(const std::string &basis_path)
{
  m_basis_path = basis_path;
}

void LagrangianInterpolation::SetRadius(const double &radius)
{
  m_radius = radius;
}

void LagrangianInterpolation::SetNumSeeds(const int &num_seeds)
{
  m_num_seeds = num_seeds;
}

void LagrangianInterpolation::SetInterval(const int &interval)
{
  m_interval = interval;
}

void LagrangianInterpolation::SetStartCycle(const int &start_cycle)
{
  m_start_cycle = start_cycle;
}

void LagrangianInterpolation::SetEndCycle(const int &end_cycle)
{
  m_end_cycle = end_cycle;
}

void LagrangianInterpolation::PreExecute()
{
  Filter::PreExecute();
  Filter::CheckForRequiredField(m_field_name);
}

void LagrangianInterpolation::PostExecute()
{
  Filter::PostExecute();
}

void LagrangianInterpolation::DoExecute()
{
#ifdef VTKH_PARALLEL
  vtkm::Id rank = vtkh::GetMPIRank();
  vtkm::Id num_ranks = vtkh::GetMPISize();
  bool allReceived[num_ranks] = {false};
  allReceived[rank] = true;

// PREPROCESSING 
// 1. READING THE FIRST FILE TO IDENTIFY BBOX + SEED VALIDITY, ADJACENT NODES 

  std::stringstream filename;
  filename << m_basis_path << "Lagrangian_Structured_" << rank << "_" << m_interval << ".vtk"; // Load the first basis flow file.

  vtkm::cont::DataSet dataset_info;
  vtkm::io::reader::VTKDataSetReader reader(filename.str().c_str());
  dataset_info = reader.ReadDataSet();

  vtkm::Bounds bounds = dataset_info.GetCoordinateSystem().GetBounds();
  double BB[6];
  BB[0] = bounds.X.Min;  
  BB[1] = bounds.X.Max;  
  BB[2] = bounds.Y.Min;  
  BB[3] = bounds.Y.Max;  
  BB[4] = bounds.Z.Min;  
  BB[5] = bounds.Z.Max;  
   
// 2. READING SEED FILE FROM DISK 

  vtkm::cont::ArrayHandle<vtkm::Particle> SeedParticleArray, SeedParticleOriginal;
  SeedParticleArray.Allocate(m_num_seeds);
  SeedParticleOriginal.Allocate(m_num_seeds);

  auto seed_portal = SeedParticleArray.GetPortalControl();
  vtkm::cont::ArrayHandle<vtkm::Id> SeedValidity;
  SeedValidity.Allocate(m_num_seeds);
  auto seed_validity = SeedValidity.GetPortalControl();

  std::ifstream seed_stream(m_seed_path);
  float x1, y1, z1; 
  int seed_counter = 0;
  
  while(seed_stream >> x1)
  {
    seed_stream >> y1;
    seed_stream >> z1;
    seed_portal.Set(seed_counter, vtkm::Particle(vtkm::Vec<vtkm::FloatDefault, 3>(x1,y1,z1),seed_counter)); 
    seed_counter++;
  }

  vtkm::cont::ArrayCopy(SeedParticleArray, SeedParticleOriginal);
  auto seed_original_portal = SeedParticleOriginal.GetPortalControl();
 
// Using extracted BBox information - check seed validity and identify seeds that belong to this block 

  vtkm::worklet::DispatcherMapField<worklets::SeedValidityCheck>(worklets::SeedValidityCheck(BB[0], BB[1], BB[2], BB[3], BB[4], BB[5])).Invoke(SeedParticleArray, SeedValidity);

// Sharing bbox information with all other nodes, so that adjacent nodes can be identified. 

  std::vector<int> neighbor_ranks;
  double *bbox_list = (double*)malloc(sizeof(double)*6*num_ranks);

  if(num_ranks > 1)
  {
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

// 2. START LOOP OVER INTERVALS 

  for(int cycle = m_start_cycle; cycle <= m_end_cycle; cycle += m_interval)
  {
    std::vector<double> px, py, pz; // Start location of a basis flow
    std::vector<double> dx, dy, dz; // Displacement of corresponding basis flow

    std::stringstream flowmap_name;
    flowmap_name << m_basis_path << "Lagrangian_Structured_" << rank << "_" << cycle << ".vtk"; // Load the first basis flow file.

    vtkm::cont::DataSet input_flowmap;
    vtkm::io::reader::VTKDataSetReader reader(flowmap_name.str().c_str());
    input_flowmap = reader.ReadDataSet();

    int num_pts = input_flowmap.GetNumberOfPoints();
    vtkm::cont::ArrayHandleVirtualCoordinates coordinatesystem = input_flowmap.GetCoordinateSystem().GetData();
    auto coords_portal = input_flowmap.GetCoordinateSystem().GetData().GetPortalControl();    
    
    auto valid_VAH = input_flowmap.GetField("valid").GetData();
    vtkm::cont::ArrayHandle<vtkm::Int32> valid_arrayhandle;
    valid_VAH.CopyTo(valid_arrayhandle);

    auto disp_VAH = input_flowmap.GetField("displacement").GetData();
    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64, 3>> disp_arrayhandle;
    disp_VAH.CopyTo(disp_arrayhandle);

    auto valid_portal = valid_arrayhandle.GetPortalControl();
    auto disp_portal = disp_arrayhandle.GetPortalControl();    


    // Perform reconstruction
    if(num_ranks)
    { 
      // Scan input_flowmap 
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
      
        // Loop over neighbor_ranks 
        // VTK objects to read flowmaps of neighboring ranks
      vtkm::cont::DataSet neighbor_flowmap;
      for(int n = 0; n < neighbor_ranks.size(); n++)
      {   
        std::stringstream neighbor_filename;
        neighbor_filename << m_basis_path << "Lagrangian_Structured_" << neighbor_ranks[n] << "_" << cycle << ".vtk";
      
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
  
      int num_basis = px.size();

      vtkm::worklet::DispatcherMapField<worklets::ShepardInterpolation>(worklets::ShepardInterpolation(m_radius, num_basis, px, py, pz, dx, dy, dz)).Invoke(valid_arrayhandle, coordinatesystem, disp_arrayhandle);
    } // End of reconstruction block
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Particle Advection Start 

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

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
      std::cout << "Completed particle advection" << std::endl;

    auto particles = res.Particles;
    auto particlePortal = particles.GetPortalControl();

    // Particle Advection End 

    // Write Output Start 

    int connectivity_index = 0;
    std::vector<vtkm::Id> connectivity;
    std::vector<vtkm::Vec<vtkm::Float64, 3>> pointCoordinates;
    std::vector<vtkm::UInt8> shapes;
    std::vector<vtkm::IdComponent> numIndices;
    std::vector<vtkm::Id> pathlineId;

    for(int i = 0; i < m_num_seeds; i++)
    {
      if(seed_validity.Get(i))
      {
        auto pt1 = seed_original_portal.Get(i).Pos;
        auto pt2 = particlePortal.Get(i).Pos;

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
      else
      {
        auto pt1 = seed_original_portal.Get(i).Pos;
        auto pt2 = particlePortal.Get(i).Pos;

        std::cout << "Invalid " << i << " " << pt1[0] << " " << pt1[1] << " " << pt1[2] << " " << pt2[0] << " " << pt2[1] << " " << pt2[2] << std::endl;
      }
    }

    vtkm::cont::DataSetBuilderExplicit DSB_Explicit;
    vtkm::cont::DataSet pathlines = DSB_Explicit.Create(pointCoordinates, shapes, numIndices, connectivity);
    vtkm::cont::DataSetFieldAdd dsfa;
    dsfa.AddCellField(pathlines, "ID", pathlineId);

    std::stringstream outputfile;
    outputfile << m_output_path << rank << "_" << cycle << ".vtk";

    vtkm::io::writer::VTKDataSetWriter writer(outputfile.str().c_str());
    writer.WriteDataSet(pathlines);
    // Write Output End

    // Exchange information - Update Seed Information Start 

    std::vector<int> outgoing_id;
    std::vector<int> outgoing_dest;

    for(int i = 0; i < m_num_seeds; i++)
    { 
      if(seed_validity.Get(i))
      { 
        auto pt = seed_portal.Get(i).Pos; 
        if(!BoundsCheck(pt[0], pt[1], pt[2], BB))
        { 
          seed_validity.Set(i, 0);
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

    if(num_ranks > 1)
    {
      MPI_Barrier(MPI_COMM_WORLD);

      int bufsize = 10*(m_num_seeds*5 + (MPI_BSEND_OVERHEAD * num_ranks)); // message packet size 4 + message packet size 1 for empty sends.
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
              auto pt = seed_portal.Get(outgoing_id[j]).Pos;
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
            seed_validity.Set(index, 1);
            seed_portal.Set(index, vtkm::Particle(vtkm::Vec<vtkm::FloatDefault, 3>
            (recvbuff[i*4+1], recvbuff[i*4+2], recvbuff[i*4+3]), index));
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
  // Exchange information - Update Seed Information End 

  } // Loop over intervals

  this->m_output = new DataSet();
  vtkm::Id domain_id;
  vtkm::cont::DataSet dom;
  this->m_input->GetDomain(0, dom, domain_id);
  m_output->AddDomain(dom, domain_id);
  
#endif
}

VTKM_EXEC
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
