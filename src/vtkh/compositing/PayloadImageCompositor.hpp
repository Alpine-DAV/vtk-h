#ifndef VTKH_DIY_PAYLOAD_IMAGE_COMPOSITOR_HPP
#define VTKH_DIY_PAYLOAD_IMAGE_COMPOSITOR_HPP

#include <vtkh/compositing/PayloadImage.hpp>
#include <algorithm>

#include<vtkh/vtkh_exports.h>

namespace vtkh
{

class VTKH_API PayloadImageCompositor
{
public:

void ZBufferComposite(vtkh::PayloadImage &front, const vtkh::PayloadImage &image)
{
  assert(front.m_depths.size() == front.m_payloads.size() / front.m_payload_bytes);
  assert(front.m_bounds.X.Min == image.m_bounds.X.Min);
  assert(front.m_bounds.Y.Min == image.m_bounds.Y.Min);
  assert(front.m_bounds.X.Max == image.m_bounds.X.Max);
  assert(front.m_bounds.Y.Max == image.m_bounds.Y.Max);

  const int size = static_cast<int>(front.m_depths.size());

#ifdef vtkh_USE_OPENMP
  #pragma omp parallel for
#endif
  for(int i = 0; i < size; ++i)
  {
    const float depth = image.m_depths[i];
    if(depth > 1.f  || front.m_depths[i] < depth)
    {
      continue;
    }
    const int offset = i * 4;
    front.m_depths[i] = depth;
    const size_t p_offset = i * front.m_payload_bytes;
    std::copy(&image.m_payloads[p_offset],
              &image.m_payloads[p_offset] + front.m_payload_bytes,
              &front.m_payloads[p_offset]);
  }
}


};

} // namespace vtkh
#endif
