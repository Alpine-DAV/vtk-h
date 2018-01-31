#ifndef VTKH_SCENE_HPP
#define VTKH_SCENE_HPP

#include <vector>
#include <list>
#include <vtkh/rendering/Render.hpp>
#include <vtkh/rendering/Renderer.hpp>

namespace vtkh 
{

class Scene 
{
private:
  std::list<vtkh::Renderer*>   m_renderers; 
  std::vector<vtkh::Render>    m_renders;
  bool                         m_has_volume;
public:
 Scene();
 ~Scene();

  void AddRender(vtkh::Render &render);
  void AddRenderer(vtkh::Renderer *render);
  void Render();
  void Save();
protected:
  bool IsMesh(vtkh::Renderer *renderer);
  bool IsVolume(vtkh::Renderer *renderer);
}; // class scene

} //namespace  vtkh
#endif
