#include <iostream>
#include <vtkh/filters/Streamline.hpp>
#include <vtkm/filter/Streamline.h>
#include <vtkh/vtkh.hpp>
#include <vtkh/Error.hpp>

#if VTKH_PARALLEL
#include <vtkm/thirdparty/diy/diy.h>
#include <vtkm/thirdparty/diy/mpi-cast.h>
#include <mpi.h>
#endif

namespace vtkh
{

Streamline::Streamline()
{
}

Streamline::~Streamline()
{

}

void Streamline::PreExecute()
{
  Filter::PreExecute();
  Filter::CheckForRequiredField(m_field_name);
}

void Streamline::PostExecute()
{
  Filter::PostExecute();
}

void Streamline::DoExecute()
{
  this->m_output = new DataSet();

#ifndef VTKH_BYPASS_VTKM_BIH

#ifdef VTKH_PARALLEL
  // Setup VTK-h and VTK-m comm.
  MPI_Comm mpi_comm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
  vtkm::cont::EnvironmentTracker::SetCommunicator(vtkmdiy::mpi::communicator(vtkmdiy::mpi::make_DIY_MPI_Comm(mpi_comm)));
#endif

  const int num_domains = this->m_input->GetNumberOfDomains();

  vtkm::cont::PartitionedDataSet inputs;
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    if(dom.HasField(m_field_name))
    {
      using vectorField_d = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64, 3>>;
      using vectorField_f = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 3>>;
      auto field = dom.GetField(m_field_name).GetData();
      if(!field.IsType<vectorField_d>() && !field.IsType<vectorField_f>())
      {
        throw Error("Vector field type does not match <vtkm::Vec<vtkm::Float32,3>> or <vtkm::Vec<vtkm::Float64,3>>");
      }
    }
    else
    {
      throw Error("Domain does not contain specified vector field for Streamline analysis.");
    }

    inputs.AppendPartition(dom);
  }

  vtkm::filter::Streamline streamlineFilter;
  auto seedsAH = vtkm::cont::make_ArrayHandle(m_seeds, vtkm::CopyFlag::Off);

  streamlineFilter.SetStepSize(m_step_size);
  streamlineFilter.SetActiveField(m_field_name);
  streamlineFilter.SetSeeds(seedsAH);
  streamlineFilter.SetNumberOfSteps(m_num_steps);
  auto out = streamlineFilter.Execute(inputs);

  for (vtkm::Id i = 0; i < out.GetNumberOfPartitions(); i++)
  {
    this->m_output->AddDomain(out.GetPartition(i), i);
  }
#endif
}

} //  namespace vtkh
