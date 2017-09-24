#ifndef VTKH_COMPOSITOR_BASE_HPP
#define VTKH_COMPOSITOR_BASE_HPP

#include <sstream>
#include <rendering/vtkh_image.hpp>

namespace vtkh 
{

class Compositor
{
public:
    enum CompositeMode { 
                         Z_BUFFER_SURFACE, // zbuffer composite no transparency
                         Z_BUFFER_BLEND,   // zbuffer composite with transparency 
                         VIS_ORDER_BLEND   // blend images in a specific order 
                       };
    Compositor();

    virtual ~Compositor();

    void SetCompositeMode(CompositeMode composite_mode);

    void ClearImages();
    
    void AddImage(const unsigned char *color_buffer,
                  const float *        depth_buffer,
                  const int            width,
                  const int            height);

    void AddImage(const float *color_buffer,
                  const float *depth_buffer,
                  const int    width,
                  const int    height);
    
    void AddImage(const unsigned char *color_buffer,
                  const int            width,
                  const int            height,
                  const int            vis_order);

    void AddImage(const float *color_buffer,
                  const int    width,
                  const int    height,
                  const int    vis_order);
    
    Image Composite();

    void SetBackgroundColor(float *background_color);

    virtual void         Cleanup();
    
    std::string          GetLogString(); 

    unsigned char * ConvertBuffer(const float *buffer, const int size)
    {
        unsigned char *ubytes = new unsigned char[size];

#ifdef VTKH_USE_OPENMP
        #pragma omp parallel for
#endif
        for(int i = 0; i < size; ++i)
        {
            ubytes[i] = static_cast<unsigned char>(buffer[i] * 255.f);
        }

        return ubytes;
    }

protected:
    virtual void CompositeZBufferSurface();
    virtual void CompositeZBufferBlend();
    virtual void CompositeVisOrder();

    std::stringstream   m_log_stream;    
    CompositeMode       m_composite_mode;
    std::vector<Image>  m_images;
    float               m_background_color[4];
};

};

#endif


