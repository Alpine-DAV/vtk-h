#ifndef VTKH_PERMUTATION_REMOVAL_HPP
#define VTKH_PERMUTATION_REMOVAL_HPP


#include <vtkm/cont/DataSet.h>

namespace vtkh {

//
// Theshold outputs CellSetPermutations which cannot
// be consumed by anything else in vtkm, so we need
// to explicitly do a deep copy and make the cell set
// explicit
//
void StripPermutation(vtkm::cont::DataSet &data_set);

}//namespace vtkh
#endif
