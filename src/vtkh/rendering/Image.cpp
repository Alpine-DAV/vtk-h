// See License.txt

#include "Image.hpp"
#include <vtkh/utils/PNGEncoder.hpp>

#include <iostream>
#include <fstream>

namespace vtkh
{

void Image::Save(const std::string &name, bool asPNG)
{
    const int width = m_bounds.X.Max - m_bounds.X.Min + 1;
    const int height = m_bounds.Y.Max - m_bounds.Y.Min + 1;

    if (asPNG)
    {
        PNGEncoder encoder;
        encoder.Encode(&m_pixels[0], width, height);
        encoder.Save(name);
    }
    else
    {
        std::ofstream outfile(name, std::ios::binary | std::ios_base::out | std::ios::trunc);
        outfile.write((char*)m_pixels.data(), width * height * 4);
        outfile.close();
    }
}

} // namespace vtkh
