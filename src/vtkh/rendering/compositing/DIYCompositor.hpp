#ifndef VTKH_DIY_COMPOSITOR_HPP
#define VTKH_DIY_COMPOSITOR_HPP

#include <vtkh/rendering/Image.hpp>
#include <vtkh/rendering/compositing/Compositor.hpp>
#include <diy/mpi.hpp>
#include <iostream>

namespace vtkh 
{

class DIYCompositor : public Compositor
{
public:
     DIYCompositor();
    ~DIYCompositor();
    
    void Cleanup() override;
    
private:
    virtual void CompositeZBufferSurface() override;
    virtual void CompositeZBufferBlend() override;
    virtual void CompositeVisOrder() override;
    diy::mpi::communicator   m_diy_comm;
    int                      m_rank;
};

}; // namespace vtkh

#endif


