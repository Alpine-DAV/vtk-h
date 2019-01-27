#include "DIYCompositor.hpp"

#include <vtkh/vtkh.hpp>
#include <vtkh/rendering/compositing/DirectSendCompositor.hpp>
#include <vtkh/rendering/compositing/RadixKCompositor.hpp>
#include <diy/mpi.hpp>

#include <assert.h>
#include <limits>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif

namespace vtkh
{
DIYCompositor::DIYCompositor()
: m_rank(0)
{
    m_diy_comm = diy::mpi::communicator(MPI_Comm_f2c(GetMPICommHandle()));
    m_rank = m_diy_comm.rank();
}

DIYCompositor::~DIYCompositor()
{
}

void
DIYCompositor::CompositeZBufferSurface()
{
  assert(m_images.size() == 1);
  RadixKCompositor compositor;

  compositor.CompositeSurface(m_diy_comm, this->m_images[0]);
  m_log_stream<<compositor.GetTimingString();

}

void
DIYCompositor::CompositeZBufferBlend()
{
  assert("this is not implemented yet" == "error");
}

void
DIYCompositor::CompositeVisOrder()
{
  assert(m_images.size() != 0);
  DirectSendCompositor compositor;
  compositor.CompositeVolume(m_diy_comm, this->m_images);
}

void
DIYCompositor::Cleanup()
{

}

}; //namespace vtkh



