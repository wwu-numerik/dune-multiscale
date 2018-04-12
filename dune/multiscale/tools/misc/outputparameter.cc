#include <config.h>

#include "outputparameter.hh"

#include <dune/multiscale/problems/selector.hh>
#include <dune/xt/common/filesystem.hh>
#include <memory>

#include "dune/multiscale/problems/base.hh"

namespace Dune {
namespace Multiscale {

// OutputParameters::OutputParameters(const std::string path
//                                  = problem.config().get("global.datadir", "data")
OutputParameters::OutputParameters(const XT::Common::Configuration& config)
  : my_prefix_("")
  , my_path_(config.get("global.datadir", "data"))
{
  Dune::XT::Common::test_create_directory(my_path_);
}

void OutputParameters::set_prefix(std::string my_prefix)
{
  my_prefix_ = my_prefix;
}

//! base of file name for data file
std::string OutputParameters::prefix() const
{
  if (my_prefix_.empty())
    return "";
  return my_prefix_ + "_";
}

//! path where the data is stored
std::string OutputParameters::path() const
{
  return my_path_;
}

int OutputParameters::outputformat() const
{
  // return 0; // GRAPE (lossless format)
  return 1; // VTK
  // return 2; // VTK vertex data
  // return 3; // gnuplot
}

bool OutputParameters::separateRankPath() const
{
  return false;
}

std::string OutputParameters::macroGridName(const int /*dim*/) const
{
  DUNE_THROW(NotImplemented, "to be removed");
  return "";
}

std::string OutputParameters::fullpath(const std::string name)
{
  std::string ret = (boost::format("%s/%s%s") % path() % prefix() % name).str();
  for (auto&& ch : ret) {
    if (ch == ' ')
      ch = '_';
  }
  return ret;
}

} // namespace Multiscale
} // namespace Dune
