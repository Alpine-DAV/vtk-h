#include <vtkh/filters/HistSampling.hpp>
#include <vtkh/filters/GhostStripper.hpp>
#include <vtkm/filter/internal/CreateResult.h>
#include <vtkm/cont/ArrayHandleTransform.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include<vtkm/io/reader/VTKDataSetReader.h>
//#include <DataSet.h>
#include <vtkm/cont/BoundsCompute.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/filter/FilterDataSetWithField.h>
#include <vtkm/filter/internal/CreateResult.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/connectivities/ImageConnectivity.h>
#include <vtkm/worklet/FieldStatistics.h>
#include<vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/cont/ArrayHandleTransform.h>
//#include <DataSet.h>
#include <vtkm/worklet/FieldHistogram.h>
//#include <vtkm/filter/ImageConnectivity.h>
#include <vtkm/filter/internal/CreateResult.h>
#include <vtkm/cont/ArrayHandleTransform.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/cont/ArrayHandleExtractComponent.h>
#include <iostream>
#include <algorithm>
#include <vtkm/worklet/WorkletMapField.h>
#ifdef VTKH_PARALLEL
  #include <mpi.h>
#endif
#include <iostream>

namespace vtkh
{

HistSampling::HistSampling()
{

}

HistSampling::~HistSampling()
{

}

void
HistSampling::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void HistSampling::PreExecute()
{
  Filter::PreExecute();
}

void HistSampling::PostExecute()
{
  Filter::PostExecute();
}

std::string
HistSampling::GetField() const
{
  return m_field_name;
}

void PrintHistogram(vtkm::cont::ArrayHandle<vtkm::Id> bins,
                    vtkm::Id numberOfBins,
                    const vtkm::Range& range,
                    vtkm::Float32 delta)
{
  vtkm::cont::ArrayHandle<vtkm::Id>::PortalConstControl binPortal = bins.GetPortalConstControl();

  vtkm::Id sum = 0;
  for (vtkm::Id i = 0; i < numberOfBins; i++)
  {
    vtkm::Float64 lo = range.Min + (static_cast<vtkm::Float64>(i) * delta);
    vtkm::Float64 hi = lo + delta;
    sum += binPortal.Get(i);
    std::cout << "  BIN[" << i << "] Range[" << lo << ", " << hi << "] = " << binPortal.Get(i)
              << std::endl;
  }
  std::cout<<"total points:"<<sum<<std::endl;
}


struct LookupWorklet : public vtkm::worklet::WorkletMapField
{
    using ControlSignature = void(FieldIn, FieldOut, WholeArrayIn);
    using ExecutionSignature = _2(_1, _3);

    template <typename TablePortal>
    VTKM_EXEC vtkm::UInt8 operator()(vtkm:: Pair<vtkm::Id,vtkm::Float32> x, TablePortal table) const {
        //return table.Get(i*i % 7);
	return x.second<table.Get(x.first);
    }
};


void PrintStatInfo(vtkm::worklet::FieldStatistics<vtkm::Float64>::StatInfo statinfo)
{
  std::cout << "   Median " << statinfo.median << std::endl;
  std::cout << "   Minimum " << statinfo.minimum << std::endl;
  std::cout << "   Maximum " << statinfo.maximum << std::endl;
  std::cout << "   Mean " << statinfo.mean << std::endl;
  std::cout << "   Variance " << statinfo.variance << std::endl;
  std::cout << "   Standard Deviation " << statinfo.stddev << std::endl;
  std::cout << "   Skewness " << statinfo.skewness << std::endl;
  std::cout << "   Kurtosis " << statinfo.kurtosis << std::endl;
  std::cout << "   Raw Moment 1-4 [ ";
  for (vtkm::Id i = 0; i < 4; i++)
    std::cout << statinfo.rawMoment[i] << " ";
  std::cout << "]" << std::endl;
  std::cout << "   Central Moment 1-4 [ ";
  for (vtkm::Id i = 0; i < 4; i++)
    std::cout << statinfo.centralMoment[i] << " ";
  std::cout << "]" << std::endl;
}

void HistSampling::DoExecute()
{
  #ifdef VTKH_PARALLEL
  this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();

  vtkh::GhostStripper stripper;

  stripper.SetInput(this->m_input);
  stripper.SetField("ascent_ghosts");
  stripper.SetMinValue(0);
  stripper.SetMaxValue(0);
  // stripper.AddMapField("point_data");

  stripper.Update();

  vtkh::DataSet *stripped_output = stripper.GetOutput();


  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::Id numberOfBins = 10;
    vtkm::Range range;
    vtkm::Float64 delta;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    stripped_output->GetDomain(i, dom, domain_id);
    // insert interesting stuff

    std::cout<<"Does it have field"<<m_field_name<<std::endl;

    if(!dom.HasField(m_field_name))
    {
      std::cout<<"Does not have field"<<m_field_name<<std::endl;
      continue;
    }


    // try to strip of ghost layers


    // std::cout << "Here:" << std::endl;
    vtkm::worklet::FieldStatistics<vtkm::Float64>::StatInfo statinfo;
    // std::cout << "Here 1:" << std::endl;
    vtkm::cont::ArrayHandle<vtkm::Float64> data;
    // std::cout << "Here 2:" << std::endl;
    dom.GetField(m_field_name).GetData().CopyTo(data);
    // std::cout << "Here 3:" << std::endl;


    vtkm::worklet::FieldStatistics<vtkm::Float64>().Run(data, statinfo);

    std::cout << "Statistics for CELL data:" << std::endl;
    PrintStatInfo(statinfo);

    vtkm::Float64 local_min = statinfo.minimum;
    vtkm::Float64 global_min = 999999999;
    vtkm::Float64 local_max = statinfo.maximum;
    vtkm::Float64 global_max = -999999999;


    // std::cout << "Here 4" << std::endl;
    MPI_Comm mpi_comm = MPI_Comm_f2c(vtkh::GetMPICommHandle());

    MPI_Allreduce(&local_min,&global_min,1, MPI_DOUBLE,MPI_MIN,mpi_comm);
    MPI_Allreduce(&local_max,&global_max,1, MPI_DOUBLE,MPI_MAX,mpi_comm);

    // vtkm::cont::ArrayHandle<vtkm::Id> binArray;
    // binArray.Allocate(numberOfBins);
    vtkm::cont::ArrayHandle<vtkm::Id> counts;
    counts.Allocate(numberOfBins);
    vtkm::worklet::FieldHistogram worklet;
    worklet.Run(data,numberOfBins,global_min,global_max,delta,counts);

    range.Min = global_min;
    range.Max = global_max;

    // range.Min = local_min;
    // range.Max = local_max;

    // std::cout << "Here 5" << "delta"<< delta<< std::endl;

    PrintHistogram(counts, numberOfBins, range, delta);
    // std::cout << "Here 6" << std::endl;

    // copy to a contig memory
    vtkm::cont::ArrayHandle<vtkm::Id>::PortalConstControl portal = counts.GetPortalConstControl();

    std::vector <vtkm::Id > myContainer(static_cast <std::size_t >( portal.GetNumberOfValues()));

    std::copy(vtkm::cont:: ArrayPortalToIteratorBegin(portal), vtkm::cont:: ArrayPortalToIteratorEnd(portal), myContainer.begin ());

    std::vector <vtkm::Id > myGlobalContainer(static_cast <std::size_t >( portal.GetNumberOfValues()));

    // std::cout << "Here 7" << std::endl;

    MPI_Allreduce(&myContainer[0],&myGlobalContainer[0],numberOfBins, MPI_INT,MPI_SUM,mpi_comm);

    MPI_Barrier(mpi_comm);

    // std::cout << "Here 8" << std::endl;

    std::copy(myGlobalContainer.begin(), myGlobalContainer.end(), std::ostream_iterator<vtkm::Float64>(std::cout, " "));

    vtkm::Id glob_sum = 0;
    for (vtkm::Id i = 0; i < numberOfBins; i++)
    {
      glob_sum += myGlobalContainer[i];
    }
    std::cout<<"total points:"<<glob_sum<<std::endl;

    vtkm::cont:: ArrayHandle <vtkm::Id > globCounts = vtkm::cont:: make_ArrayHandle(myGlobalContainer);

    // MPI_Barrier(mpi_comm);

    // MPI_Finalize();

    // start doing sampling

    float sample_percent = 0.1; // To-Do make it an argument

    vtkm::cont:: ArrayHandle <vtkm::Float64> output;

    //auto field = dom.GetField(m_field_name);

    // vtkm::cont:: ArrayHandleExtractComponent <vtkm::cont:: ArrayHandle <vtkm::Vec <vtkm::Float64 , 3>> > field(field_orig , 0);

    vtkm::Int32 tot_points = data.GetNumberOfValues();

    vtkm::cont:: ArrayHandleIndex indexArray (numberOfBins);
    vtkm::cont::ArrayHandle<vtkm::Id> indices;
    vtkm::cont:: Algorithm ::Copy(indexArray,indices);

    vtkm::cont:: ArrayHandleZip <vtkm::cont:: ArrayHandle <vtkm::Id >, vtkm::cont:: ArrayHandle <vtkm::Id >> zipArray(globCounts , indices );  // change for global bin Array

    vtkm::cont:: Algorithm ::Sort(zipArray );
    // vtkm::cont::ArrayHandle<vtkm::Pair< vtkm::Id , vtkm::Id>> dataVals;
    // zipArray.CopyTo(dataVals);

    vtkm::cont:: ArrayHandleZip <vtkm::cont:: ArrayHandle <vtkm::Id >, vtkm::cont:: ArrayHandle <vtkm::Id >>::PortalConstControl binPortal = zipArray.GetPortalConstControl();

    for (int i = 0; i < numberOfBins; ++i)
    {
      std::cout<<binPortal.Get(i)<<std::endl;
    }

    // create a float array now to hold the fractional counts for computing acceptance histogram
    vtkm::cont::ArrayHandle<vtkm::Float32> targetCounts;

    vtkm::Float32 remainingSamples = sample_percent*tot_points;
    //cout<<"sample_percent "<<sample_percent<<" tot_points "<<tot_points<<endl;
    vtkm::Float32 remainingBins = numberOfBins;
    std::vector<vtkm::Float32> targetSamples;

    for (int i = 0; i < numberOfBins; ++i)
    {
      vtkm::Float32 targetNeededSamples = remainingSamples/(1.0*remainingBins);
      vtkm::Float32 curCount = (vtkm::Float32)binPortal.Get(i).first;
      vtkm::Float32 samplesTaken;
      // cout<<"curCount "<<curCount<<" targetNeededSamples "<<targetNeededSamples<<endl;

      if (curCount<targetNeededSamples)
      {
        samplesTaken = curCount;
      }
      else // for speed up, this else loop can be used to set the rest of the samples
      {
        samplesTaken = targetNeededSamples;
      }
      targetSamples.push_back(samplesTaken);
      remainingBins = remainingBins-1;
      remainingSamples = remainingSamples - samplesTaken;
    }
    std::vector<vtkm::Float32> acceptanceProbsVec(numberOfBins,-1.0);
    vtkm::Float32 sum=0.0;
    int counter=0;
    for(vtkm::Float32 n : targetSamples)
    {
      //acceptanceProbsVec.push_back(n/binPortal.Get(counter).first);
      acceptanceProbsVec[binPortal.Get(counter).second] = n/binPortal.Get(counter).first;
      std::cout<<"samples:"<<n<<" prob:"<<n/binPortal.Get(counter).first<<std::endl;
      sum+=n;
      counter++;

    }
    counter=0;
    for(vtkm::Float32 n : acceptanceProbsVec)
    {
      std::cout<<"prob:"<<n<<std::endl;
    }
    std::cout<<"Samples requested:"<<sample_percent*tot_points<<std::endl;
    std::cout<<"Samples expected:"<<sum<<std::endl;

    // use the acceptance probabilities to create a stencil buffer
    struct convertToBinId
    {
      VTKM_CONT convertToBinId(vtkm::Id nbins, vtkm::Range range):numberOfBins(nbins),rangeOfValues(range){}
      VTKM_EXEC_CONT  vtkm::Id operator() ( vtkm:: Float32 x)
      const
      {
        vtkm::Id id = (vtkm::Id)(numberOfBins* (x - rangeOfValues.Min) / ( rangeOfValues.Max -rangeOfValues.Min));
        if (id==numberOfBins)
        {
          id--;
        }
        return id;
      }

      // vtkm::Float32 sample_percent;
      vtkm::Id numberOfBins;
      vtkm::Range rangeOfValues;
    };

    // vtkm::cont:: ArrayHandleTransform<vtkm::cont:: ArrayHandle <vtkm::Float32 >, convertToBinId> (dataArray, convertToBinId(numberOfBins, range));

    vtkm::cont::ArrayHandle<vtkm::Id> binIndex;
      binIndex.Allocate(tot_points);

      // Worklet to set the bin number for each data value
      vtkm::worklet::FieldHistogram::SetHistogramBin<vtkm::Float64> binWorklet(numberOfBins, range.Min, delta);

      vtkm::worklet::DispatcherMapField <vtkm::worklet::FieldHistogram::SetHistogramBin<vtkm::Float64>>
        setHistogramBinDispatcher(binWorklet);

      setHistogramBinDispatcher.Invoke(data, binIndex);

    /*  struct equalsK1
    {
      equalsK1(vtkm::cont:: ArrayHandle <vtkm::Float32 > sample_probs):sample_probabilites(sample_probs){}
      equalsK1(){};
      VTKM_EXEC_CONT  bool operator() ( vtkm:: Pair<vtkm::Id,vtkm::Float32> x)
      const
      {
        vtkm::cont::ArrayHandle<vtkm::Float32>::PortalConstControl probPortal = sample_probabilites.GetPortalConstControl();
        //cout<<x.first<<" "<<x.second<<" "<<probPortal.Get(x.first)<<endl;
        return x.second<probPortal.Get(x.first);
      }


      vtkm::cont:: ArrayHandle <vtkm::Float32 > sample_probabilites;
    };

    struct equalsK2
    {
      equalsK2(vtkm::cont:: ArrayHandle <vtkm::Float32 > sample_probs):sample_probabilites(sample_probs){}
      VTKM_EXEC_CONT  bool operator() ( vtkm:: Pair<vtkm::Id,vtkm::Float32> x)
      const
      {
        vtkm::cont::ArrayHandle<vtkm::Float32>::PortalConstControl probPortal = sample_probabilites.GetPortalConstControl();
        //cout<<x.first<<" "<<x.second<<" "<<probPortal.Get(x.first)<<endl;
        return x.second<probPortal.Get(x.first);
      }

      vtkm::cont:: ArrayHandle <vtkm::Float32 > sample_probabilites;
    };*/

    // struct equalsK
    // {
    //  equalsK(vtkm::Float32 samp_percent):sample_percent(samp_percent){}
    //  VTKM_EXEC_CONT  bool operator() ( vtkm:: Float32 x) const { return x<sample_percent; }

    //  vtkm::Float32 sample_percent;
    // };

    struct randArraygen
    {
      VTKM_EXEC_CONT  vtkm::Float32 operator() ( vtkm::UInt32 seed)
    const  // RAND_MAX assumed to be 32767
      {
          seed = seed * 1103515245 + 12345;
          return ((unsigned int)(seed/65536) % 32768)/32767.0;
      }
    };

    struct randArraygen2
    {
      VTKM_EXEC_CONT  vtkm::Float32 operator() ( vtkm::UInt32 seed)
        const
          {
            vtkm::Float32 x = 0.0f;
            vtkm::Float32 xadd = 1.0f;
            vtkm::UInt32 b2 = 1 + seed;
            while (b2 != 0)
            {
              xadd *= 0.5f;
              if ((b2 & 1) != 0)
                x += xadd;
              b2 >>= 1;
            }
            return x;
          }
    };

    // vtkm::cont::ArrayHandleImplicit<randArraygen> randArray (randArraygen(),tot_points);
    //vtkm::cont::ArrayHandleImplicit<randArraygen2> randArray (randArraygen2(),tot_points);

    vtkm::cont::ArrayHandle<vtkm::Float32> randArray;

   randArray.Allocate(tot_points);

   auto randPortal = randArray.GetPortalControl();

   for(int i=0; i<tot_points; i++)
   {
	randArraygen2 randgen;
	randPortal.Set(i, randgen(i));
   }

    auto stencil = vtkm::cont:: make_ArrayHandleZip(binIndex , randArray);

    vtkm::cont:: ArrayHandle <vtkm:: Pair<vtkm::Id,vtkm::Float32>> stencil2;

    vtkm::cont:: Algorithm ::Copy(stencil , stencil2 );


    //convert vec to ArrayHandle from acceptanceProbsVec
    vtkm::cont:: ArrayHandle <vtkm::Float32 > probArray = vtkm::cont:: make_ArrayHandle(acceptanceProbsVec);

    vtkm::cont:: ArrayHandle <vtkm::UInt8> stencilBool;

    vtkm::worklet::DispatcherMapField<LookupWorklet> dispatcher;
    dispatcher.Invoke(stencil2, stencilBool, probArray);

    //auto stencilBool = vtkm::cont:: make_ArrayHandleTransform(stencil , equalsK1(probArray));


    vtkm::cont:: Algorithm ::Copy(stencilBool , output );

    vtkm::cont:: DataSetFieldAdd dataSetFieldAdd;

    dataSetFieldAdd.AddPointField(dom , "valSampled", output );



    // end sampling

    // #endif

    // std::cout << "Here 9" << std::endl;

    // end interesting stuff

    //m_output->AddDomain(valsampled, domain_id);

    m_output->AddDomain(dom, domain_id);
  }
  #endif
}

std::string
HistSampling::GetName() const
{
  return "vtkh::HistSampling";
}

} //  namespace vtkh
