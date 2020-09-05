#include <vtkh/filters/HistSampling.hpp>
#include <vtkh/filters/GhostStripper.hpp>
#include <vtkh/filters/Threshold.hpp>
#include <vtkh/filters/Histogram.hpp>
#include <vtkh/filters/ValGradHistogram.hpp>
#include <vtkh/Error.hpp>

#include <vtkh/filters/Gradient.hpp>
#include <vtkh/filters/VectorMagnitude.hpp>


#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>

#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/worklet/FieldStatistics.h>
#include <vtkm/filter/CreateResult.h>
#include <vtkm/cont/ArrayHandleTransform.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <iostream>
#include <algorithm>
#include <vtkm/worklet/WorkletMapField.h>
#include <iostream>

namespace vtkh
{

namespace detail
{

class RandomGenerate : public vtkm::worklet::WorkletMapField
{
protected:
  vtkm::Int32 m_seed;
public:
  VTKM_CONT
  RandomGenerate(vtkm::Int32 seed)
   : m_seed(seed)
  {}

  typedef void ControlSignature(FieldOut);
  typedef void ExecutionSignature(WorkIndex, _1);

  VTKM_EXEC
  void operator()(const vtkm::Id &index, vtkm::Float32 &value) const
  {
    const vtkm::Int32 sample = static_cast<vtkm::UInt32>(m_seed + index);
    vtkm::Float32 y = 0.0f;
    vtkm::Float32 yadd = 1.0f;
    vtkm::Int32 bn = sample;
    const vtkm::Int32 base = 7;
    while (bn != 0)
    {
      yadd *= 1.0f / (vtkm::Float32)base;
      y += (vtkm::Float32)(bn % base) * yadd;
      bn /= base;
    }

    value = y;
  }
}; //class RandomGenerate


vtkm::cont::ArrayHandle<vtkm::Float32>
calculate_1d_pdf(const vtkm::Int32 tot_points,
              const vtkm::Int32 num_bins,
              const vtkm::Float32 sample_percent,
              vtkm::cont::ArrayHandle<vtkm::Id> mybins)
{
  vtkm::cont:: ArrayHandle <vtkm::Id > bins;
  vtkm::cont:: Algorithm ::Copy(mybins , bins);
  vtkm::cont::ArrayHandleIndex indexArray (num_bins);
  vtkm::cont::ArrayHandle<vtkm::Id> indices;
  vtkm::cont::Algorithm::Copy(indexArray, indices);

  vtkm::cont:: ArrayHandleZip <vtkm::cont:: ArrayHandle <vtkm::Id >,
                               vtkm::cont:: ArrayHandle <vtkm::Id >>
                                 zipArray(bins, indices );

  vtkm::cont::Algorithm::Sort(zipArray);

  auto binPortal = zipArray.ReadPortal();

  vtkm::Float32 remainingSamples = sample_percent*tot_points;

  vtkm::Float32 remainingBins = num_bins;
  std::vector<vtkm::Float32> targetSamples;

  for (int i = 0; i < num_bins; ++i)
  {
    vtkm::Float32 targetNeededSamples = remainingSamples / (1.0f*remainingBins);
    vtkm::Float32 curCount = (vtkm::Float32)binPortal.Get(i).first;
    vtkm::Float32 samplesTaken;

    if(curCount < targetNeededSamples)
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

  vtkm::cont::ArrayHandle<vtkm::Float32> acceptanceProbsVec;
  acceptanceProbsVec.Allocate(num_bins);
  auto acceptance_portal = acceptanceProbsVec.WritePortal();
  for(int i = 0; i < num_bins; ++i)
  {
    acceptance_portal.Set(i, -1.f);
  }

  vtkm::Float32 sum=0.0;
  int counter=0;
  for(vtkm::Float32 n : targetSamples)
  {
    acceptance_portal.Set(binPortal.Get(counter).second,n/binPortal.Get(counter).first);
    if (binPortal.Get(counter).first < 0.00000000000001f)
    {
    	acceptance_portal.Set(binPortal.Get(counter).second,0.0);
    }
    else
    {
      acceptance_portal.Set(binPortal.Get(counter).second,n/binPortal.Get(counter).first);
    }
    sum+=n;
    counter++;

  }
  counter = 0;

  return acceptanceProbsVec;
}

vtkm::cont::ArrayHandle<vtkm::Float32>
calculate_2d_pdf(const vtkm::Int32 tot_points,
              const vtkm::Float32 sample_percent,
              ValGradHistogram::NonSparseHistogramResult bivar_histogram)
{

  vtkm::Int32 num_bins = bivar_histogram.m_nob;

  Histogram::HistogramResult valMarginal = bivar_histogram.getValDist();

  vtkm::cont::ArrayHandleIndex indexArray(num_bins);
  vtkm::cont::ArrayHandle<vtkm::Id> indices;
  vtkm::cont::Algorithm::Copy(indexArray,indices);

  vtkm::cont::ArrayHandle<vtkm::Id> val_counts_copy;
  vtkm::cont::Algorithm::Copy(valMarginal.m_bins,val_counts_copy);

  vtkm::cont::ArrayHandleZip <vtkm::cont:: ArrayHandle <vtkm::Id >, 
                              vtkm::cont:: ArrayHandle <vtkm::Id >> 
                                zipArray(val_counts_copy , indices );

  //sorts the histogram bins based on their counts
  vtkm::cont:: Algorithm ::Sort(zipArray );
  vtkm::cont::ArrayHandleZip <vtkm::cont::ArrayHandle <vtkm::Id >, 
                              vtkm::cont:: ArrayHandle <vtkm::Id >>
                                ::PortalConstControl binPortal = zipArray.GetPortalConstControl();

  

  vtkm::Float32 remainingSamples = sample_percent*tot_points;
  vtkm::Float32 remainingBins = num_bins;
  std::vector<vtkm::Float32> targetSamples(num_bins, 0.0);

  for (int i = 0; i < num_bins; ++i)
  {
    vtkm::Float32 targetNeededSamples = remainingSamples/(1.0*remainingBins);
    vtkm::Float32 curCount = (vtkm::Float32)binPortal.Get(i).first;
    vtkm::Float32 samplesTaken;
    
    if (curCount<targetNeededSamples)
    {
      samplesTaken = curCount;
    }
    else
    {
      samplesTaken = targetNeededSamples;
    }
    targetSamples[binPortal.Get(i).second] = samplesTaken;
    remainingBins = remainingBins-1;
    remainingSamples = remainingSamples - samplesTaken; 
  }


  // Start calculating the 2d acceptance probability histogram (flattened)
  vtkm::cont::ArrayHandle<vtkm::Float32> acceptanceProbsVec;
  acceptanceProbsVec.Allocate(num_bins*num_bins);
  auto acceptance_portal = acceptanceProbsVec.WritePortal();
  for(int i = 0; i < num_bins*num_bins; ++i)
  {
    acceptance_portal.Set(i, -1.f);
  }
  

  // for each bin in the val histogram, identify how many samples to pick from each 
  // bin of the corresponding marginal gradient histogram bin

  vtkm::cont::ArrayHandle<vtkm::Id> biVarFreq;
  vtkm::cont::Algorithm::Copy(bivar_histogram.m_freqs,biVarFreq);
  auto biVarFreqPortal = biVarFreq.ReadPortal();

  for(int i=0; i<num_bins; i++)
  {
    vtkm::Float32 current_bin_target_samples = targetSamples[i];

    for(int j=num_bins-1; j>=0 ; j--)
    {
      vtkm::Id curGradFreq = biVarFreqPortal.Get(i*num_bins + j);
      if(current_bin_target_samples <= curGradFreq)
      {
        acceptance_portal.Set(i*num_bins + j, current_bin_target_samples);
        current_bin_target_samples = 0;
      }
      else{
        acceptance_portal.Set(i*num_bins + j, curGradFreq);
        current_bin_target_samples -= curGradFreq;
      }
      if(curGradFreq > 0)
      {
        vtkm::Float32 temp_val = (acceptance_portal.Get(i*num_bins + j)*1.0)/curGradFreq;
        acceptance_portal.Set(i*num_bins + j, temp_val);
      }
      else
        acceptance_portal.Set(i*num_bins + j, 0.0);
    }
  }

  return acceptanceProbsVec;
}


}

HistSampling::HistSampling()
  : m_sample_percent(0.1f),
    m_num_bins(128),
    m_use_gradient(false)
{

}

HistSampling::~HistSampling()
{

}

void
HistSampling::SetSamplingPercent(const float percent)
{
  if(percent <= 0.f || percent > 1.f)
  {
    throw Error("HistSampling: sampling percent must be in the range (0,1]");
  }
  m_sample_percent = percent;
}

void
HistSampling::SetNumBins(const int num_bins)
{
  if(num_bins <= 0)
  {
    throw Error("HistSampling: num_bins must be positive");
  }
  m_num_bins = num_bins;
}

void
HistSampling::SetGradientSampling(const bool use_gradient)
{
  m_use_gradient = use_gradient;
}

void
HistSampling::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void
HistSampling::SetGhostField(const std::string &field_name)
{
  m_ghost_field = field_name;
}

void HistSampling::PreExecute()
{
  Filter::PreExecute();
  Filter::CheckForRequiredField(m_field_name);
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

struct Lookup_1d_Worklet : public vtkm::worklet::WorkletMapField
{
protected:
  vtkm::Id m_num_bins;
  vtkm::Float64 m_min;
  vtkm::Float64 m_bin_delta;
public:
  Lookup_1d_Worklet(const vtkm::Id num_bins,
                const vtkm::Float64 min_value,
                const vtkm::Float64 bin_delta)
    : m_num_bins(num_bins),
      m_min(min_value),
      m_bin_delta(bin_delta)
  {}

  using ControlSignature = void(FieldIn, FieldOut, WholeArrayIn, FieldIn);
  using ExecutionSignature = _2(_1, _3, _4);

  template <typename TablePortal>
  VTKM_EXEC vtkm::UInt8 operator()(const vtkm::Float64 &field_value,
                                   TablePortal table,
                                   const vtkm::Float32 &random) const
  {
    vtkm::Id bin = static_cast<vtkm::Id>((field_value - m_min) / m_bin_delta);
    if(bin < 0)
    {
      bin = 0;
    }
    if(bin >= m_num_bins)
    {
      bin = m_num_bins - 1;
    }

    return random < table.Get(bin);
  }
};

struct Lookup_2d_Worklet : public vtkm::worklet::WorkletMapField
{
protected:
  vtkm::Id m_num_bins;
  vtkm::Float64 m_min_val;
  vtkm::Float64 m_bin_delta_val;
  vtkm::Float64 m_min_grad;
  vtkm::Float64 m_bin_delta_grad;
public:
  Lookup_2d_Worklet(const vtkm::Id num_bins,
                const vtkm::Float64 min_val,
                const vtkm::Float64 delta_val,
                const vtkm::Float64 min_grad,
                const vtkm::Float64 delta_grad)
    : m_num_bins(num_bins),
      m_min_val(min_val),
      m_bin_delta_val(delta_val),
      m_min_grad(min_grad),
      m_bin_delta_grad(delta_grad)
  {}

  using ControlSignature = void(FieldIn, FieldIn, FieldOut, WholeArrayIn, FieldIn);
  using ExecutionSignature = _3(_1, _2, _4, _5);

  template <typename TablePortal>
  VTKM_EXEC vtkm::UInt8 operator()(const vtkm::Float64 &field_val,
                                   const vtkm::Float64 &field_grad,
                                   TablePortal table,
                                   const vtkm::Float32 &random) const
  {
    vtkm::Id bin_val = static_cast<vtkm::Id>((field_val - m_min_val) / m_bin_delta_val);
    if(bin_val < 0)
    {
      bin_val = 0;
    }
    if(bin_val >= m_num_bins)
    {
      bin_val = m_num_bins - 1;
    }

    vtkm::Id bin_grad = static_cast<vtkm::Id>((field_grad - m_min_grad) / m_bin_delta_grad);
    if(bin_grad < 0)
    {
      bin_grad = 0;
    }
    if(bin_grad >= m_num_bins)
    {
      bin_grad = m_num_bins - 1;
    }

    //flattened bin index
    vtkm::Id bin = bin_val*m_num_bins + bin_grad;

    return random < table.Get(bin);
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
  vtkh::DataSet *input = this->m_input;
  bool has_ghosts = m_ghost_field != "";

  if(has_ghosts)
  {
    vtkh::GhostStripper stripper;

    stripper.SetInput(this->m_input);
    stripper.SetField(m_ghost_field);
    stripper.SetMinValue(0);
    stripper.SetMaxValue(0);
    stripper.Update();
    input = stripper.GetOutput();
  }

   
  vtkm::cont::ArrayHandle <vtkm::Float32 > probArray;
  vtkm::Range val_range, grad_range;
  vtkm::Float64 val_delta, grad_delta;
  const int num_domains = input->GetNumberOfDomains();
  vtkh::DataSet *mag_output;


  

  if(!m_use_gradient)
  {
    Histogram histogrammer;
    histogrammer.SetNumBins(m_num_bins);

    Histogram::HistogramResult histogram = histogrammer.Run(*input, m_field_name);
    // histogram.Print(std::cout);

    vtkm::Id global_num_values = histogram.totalCount();
    vtkm::cont:: ArrayHandle <vtkm::Id > globCounts = histogram.m_bins;

    val_range = histogram.m_range;
    val_delta = histogram.m_bin_delta;

    // calculate the acceptance probability histogram 
    probArray = detail::calculate_1d_pdf(global_num_values, m_num_bins, m_sample_percent, globCounts);
  }
  else
  {
    //Calculate the gradient magnitude field
    vtkh::Gradient grad;
    grad.SetInput(input);
    grad.SetField(m_field_name);

    vtkh::GradientParameters params;
    params.output_name = "grad";
    params.use_point_gradient = false;
    grad.SetParameters(params);

    grad.Update();

    vtkh::DataSet *grad_output = grad.GetOutput();

    vtkh::VectorMagnitude mag;
    mag.SetInput(grad_output);
    mag.SetField("grad");

    mag.SetResultName("mag");
    mag.Update();

    mag_output = mag.GetOutput();


    const int grad_num_domains = mag_output->GetNumberOfDomains();
    //make sure both the fields have same number of domains
    if( num_domains != grad_num_domains)
    {
      throw Error("Number of domains does not match for the scalar value field and the gradient magnitude field");
    }

    ValGradHistogram vg_histogrammer;
    vg_histogrammer.SetNumBins(m_num_bins);

    ValGradHistogram::NonSparseHistogramResult bivar_histogram = vg_histogrammer.Run(*input, m_field_name, *mag_output);
    // bivar_histogram.Print(std::cout);

    vtkm::Id global_num_values = bivar_histogram.totalCount();

    val_range = bivar_histogram.m_range_val;
    val_delta = bivar_histogram.m_bin_delta_val;
    grad_range = bivar_histogram.m_range_grad;
    grad_delta = bivar_histogram.m_bin_delta_grad;

    // calculate the acceptance probability histogram 
    probArray = detail::calculate_2d_pdf(global_num_values, m_sample_percent, bivar_histogram);
  }

  bool valid_field;
  vtkm::cont::Field::Association assoc = input->GetFieldAssociation(m_field_name, valid_field);
  
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::cont::DataSet &dom = input->GetDomain(i);

    if(!dom.HasField(m_field_name))
    {
      // We have already check to see if the field exists globally,
      // so just skip if this particular domain doesn't have the field
      continue;
    }

    vtkm::cont::ArrayHandle<vtkm::Float64> val_data;
    dom.GetField(m_field_name).GetData().CopyTo(val_data);

    vtkm::Int32 tot_points = val_data.GetNumberOfValues();

    // use the acceptance probabilities to create a stencil buffer
    vtkm::cont::ArrayHandle<vtkm::Float32> randArray;
    randArray.Allocate(tot_points);
    const vtkm::Int32 seed = 0;
    vtkm::worklet::DispatcherMapField<detail::RandomGenerate>(seed).Invoke(randArray);

    vtkm::cont::ArrayHandle <vtkm::UInt8> stencilBool;

    if(!m_use_gradient)
    {
      vtkm::worklet::DispatcherMapField<Lookup_1d_Worklet>(Lookup_1d_Worklet{m_num_bins,
                                                       val_range.Min,
                                                       val_delta}).Invoke(val_data,
                                                                          stencilBool,
                                                                          probArray,
                                                                          randArray);
    }
    else
    {
      vtkm::cont::DataSet &grad_dom = mag_output->GetDomain(i);

      if(!grad_dom.HasField("mag"))
      {
        // We have already check to see if the field exists globally,
        // so just skip if this particular domain doesn't have the field
        continue;
      }

      vtkm::cont::ArrayHandle<vtkm::Float64> grad_data;
      grad_dom.GetField("mag").GetData().CopyTo(grad_data);

      if(tot_points != grad_data.GetNumberOfValues())
      {
        //Sanity check: the number of points in the scalar value field 
        // and the gradient magnitude field must match.
        throw Error("Mismatch in the number of values for scalar values and gradient magnitude");
      }

      
      vtkm::worklet::DispatcherMapField<Lookup_2d_Worklet>(Lookup_2d_Worklet{m_num_bins,
                                                         val_range.Min,
                                                         val_delta,
                                                         grad_range.Min,
                                                         grad_delta}).Invoke(val_data,
                                                                            grad_data,
                                                                            stencilBool,
                                                                            probArray,
                                                                            randArray);
    }

    vtkm::cont::ArrayHandle <vtkm::Float32> output;
    vtkm::cont::Algorithm ::Copy(stencilBool , output);

    // // Test code: Verify the number of eventually selected samples by summing the stencilBool
    // vtkm::Float32 test_sum = 0.0;
    // for(int j=0; j<tot_points; j++)
    // {
    //   test_sum += output.ReadPortal().Get(j);
    // }
    // std::cout << "Samples taken = " << test_sum << ". Total points = " 
    //           << tot_points << ". At domain = "<< i << std::endl;

    vtkm::cont:: DataSetFieldAdd dataSetFieldAdd;

    if(assoc == vtkm::cont::Field::Association::POINTS)
    {
      dataSetFieldAdd.AddPointField(dom , "valSampled", output );
    }
    else
    {
      dataSetFieldAdd.AddCellField(dom , "valSampled", output );
    }
  }

  vtkh::Threshold thresher;
  thresher.SetInput(input);
  thresher.SetField("valSampled");

  double upper_bound = 1.;
  double lower_bound = 1.;

  thresher.SetUpperThreshold(upper_bound);
  thresher.SetLowerThreshold(lower_bound);
  thresher.Update();
  this->m_output = thresher.GetOutput();
  
  if(has_ghosts)
  {
    delete input;
  }
  if(m_use_gradient){
    delete mag_output;
  }
}

std::string
HistSampling::GetName() const
{
  return "vtkh::HistSampling";
}

} //  namespace vtkh
