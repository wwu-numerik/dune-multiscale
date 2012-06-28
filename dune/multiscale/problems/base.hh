#ifndef DUNE_MS_PROBLEMS_BASE_HH
#define DUNE_MS_PROBLEMS_BASE_HH

#include <string>
#include <dune/multiscale/problems/constants.hh>

namespace Problem {
class IModelProblemData
{
protected:
  // name of the file where data is saved
  const std::string file_name_;
  int current_number_of_cell_problem_;
  const Constants constants_;

public:
  // Constructor for ModelProblemData
  inline IModelProblemData(const Constants constants, const std::string file_name = "no_name")
    : file_name_(file_name)
      , constants_(constants)
      , current_number_of_cell_problem_(-1)
  {}

  // epsilon (the smaller epsilon, the finer the micro-structure)
  // in the periodic setting, epsilon denotes the periode of the fine-scale oscillations
  // in the non-periodic setting, can be seen as a representative size for the fine-scale behaviour
  inline double getEpsilon() const {
    return constants_.epsilon;
  }

  // epsilon (the smaller epsilon, the finer the micro-structure)
  inline double getEpsilonEstimated() const {
    return constants_.epsilon_est;
  }

  // edge length of a cell (where we solve the cell problems)
  // we need delta >= epsilon
  inline double getDelta() const {
    return constants_.delta;
    // NOTE that (delta/epsilon_est) needs to be a positive integer!
  }

  // get an information on whether we use the solutions of cell problems that are already computed and saved in a file
  // with the name 'name_'
  inline std::string getName_and_getBool(bool& use_saved) const {
    use_saved = (file_name_ != "no_name");
    return file_name_;
  }

  inline void set_current_number_of_cell_problem(int number) {
    current_number_of_cell_problem_ = number;
  }

  inline int get_current_number_of_cell_problem() {
    return current_number_of_cell_problem_;
  }

  virtual int  getRefinementLevelReferenceProblem() const = 0;
  virtual void getMacroGridFile(std::string& macroGridName) const = 0;
  virtual int  get_Number_of_Model_Problem() const = 0;
};
}

#endif // DUNE_MS_PROBLEMS_BASE_HH
