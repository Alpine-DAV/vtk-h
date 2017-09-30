// See License.txt

#include "Image.hpp"
#include <vtkh/utils/PNGEncoder.hpp>

namespace vtkh
{

void Image::Save(std::string name)
{
    PNGEncoder encoder;
    encoder.Encode(&m_pixels[0],
        m_bounds.X.Max - m_bounds.X.Min + 1,
        m_bounds.Y.Max - m_bounds.Y.Min + 1);
    encoder.Save(name);
}

} // namespace vtkh
