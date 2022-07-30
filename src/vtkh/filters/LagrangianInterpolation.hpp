#ifndef VTK_H_LAGRANGIANINTERPOLATION_HPP
#define VTK_H_LAGRANGIANINTERPOLATION_HPP

#include <vtkh/vtkh_exports.h>
#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>

namespace vtkh
{

class VTKH_API LagrangianInterpolation : public Filter
{
public:
  LagrangianInterpolation();
  virtual ~LagrangianInterpolation();
  std::string GetName() const override;
	void SetField(const std::string &field_name);
	void SetSeedPath(const std::string &seed_path);
	void SetOutputPath(const std::string &output_path);
	void SetBasisPath(const std::string &basis_path);
	void SetRadius(const double &radius);
	void SetNumSeeds(const int &num_seeds);
	void SetInterval(const int &interval);
	void SetStartCycle(const int &start_cycle);
	void SetEndCycle(const int &end_cycle);

private:
  bool AllMessagesReceived(bool *a, int num_ranks);
  bool BoundsCheck(float x, float y, float z, double *BB);
  bool BoundsCheck(float x, float y, float z, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
  std::vector<int> GetNeighborRankList(vtkm::Id num_ranks, vtkm::Id rank, double *bbox_list);



protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  std::string m_field_name;
	double m_step_size;
	int m_write_frequency;
	int m_cust_res;
	int m_x_res, m_y_res, m_z_res;
  std::string m_seed_path;
  std::string m_basis_path;
  std::string m_output_path;
  double m_radius;
  int m_num_seeds;
  int m_interval;
  int m_start_cycle;
  int m_end_cycle;


};

} //namespace vtkh
#endif
