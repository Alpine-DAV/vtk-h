//-----------------------------------------------------------------------------
///
/// file: vtkm_permutation_removal.cpp
///
//-----------------------------------------------------------------------------


#include <vtkh/utils/vtkm_permutation_removal.hpp>
#include <vtkm/cont/CellSetPermutation.h>
#include <vtkm/worklet/CellDeepCopy.h>
namespace vtkh
{

typedef vtkm::cont::CellSetPermutation<vtkm::cont::CellSetStructured<2>>
  PermStructured2d;

typedef vtkm::cont::CellSetPermutation<vtkm::cont::CellSetStructured<3>>
  PermStructured3d;

typedef vtkm::cont::CellSetPermutation<vtkm::cont::CellSetExplicit<>>
  PermExplicit;

typedef  vtkm::cont::CellSetPermutation<vtkm::cont::CellSetSingleType<>>
  PermExplicitSingle;

void StripPermutation(vtkm::cont::DataSet &data_set)
{
  vtkm::cont::DynamicCellSet cell_set = data_set.GetCellSet();
  vtkm::cont::DataSet result;
  vtkm::cont::CellSetExplicit<> explicit_cells;

  if(cell_set.IsSameType(PermStructured2d()))
  {
    PermStructured2d perm = cell_set.Cast<PermStructured2d>();
    explicit_cells = vtkm::worklet::CellDeepCopy::Run(perm);
  }
  else if(cell_set.IsSameType(PermStructured3d()))
  {
    PermStructured3d perm = cell_set.Cast<PermStructured3d>();
    explicit_cells = vtkm::worklet::CellDeepCopy::Run(perm);
  }
  else if(cell_set.IsSameType(PermExplicit()))
  {
    PermExplicit perm = cell_set.Cast<PermExplicit>();
    explicit_cells = vtkm::worklet::CellDeepCopy::Run(perm);
  }
  else if(cell_set.IsSameType(PermExplicitSingle()))
  {
    PermExplicitSingle perm = cell_set.Cast<PermExplicitSingle>();
    explicit_cells = vtkm::worklet::CellDeepCopy::Run(perm);
  }

  result.AddCellSet(explicit_cells);

  vtkm::Id num_coords = data_set.GetNumberOfCoordinateSystems();
  for(vtkm::Id i = 0; i < num_coords; ++i)
  {
    result.AddCoordinateSystem(data_set.GetCoordinateSystem(i));
  }

  vtkm::Id num_fields = data_set.GetNumberOfFields();
  for(vtkm::Id i = 0; i < num_fields; ++i)
  {
    result.AddField(data_set.GetField(i));
  }

  data_set = result;
}
} // namespace vtkh

