#ifndef VTK_H_FILTER_HPP
#define VTK_H_FILTER_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/DataSet.hpp>

namespace vtkh
{

class Filter
{
public:
  Filter();
  virtual ~Filter();
  void SetInput(DataSet *input);
  virtual std::string GetName() const = 0;

  DataSet* GetOutput();
  DataSet* Update();

  void AddMapField(const std::string &field_name);

  void ClearMapFields();

protected:
  virtual void DoExecute() = 0;
  virtual void PreExecute();

  virtual void PostExecute();

  std::vector<std::string> m_map_fields;

  DataSet *m_input;
  DataSet *m_output;

  void MapAllFields();

  void PropagateMetadata();
};

} //namespace vtkh
#endif
