// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef OUTPUTPARAMETER_HH
#define OUTPUTPARAMETER_HH

#include <string>

#include <dune/stuff/common/parameter/configcontainer.hh>

namespace Dune {
namespace Multiscale {
//! define output parameters for \ref Dune::DataOutput
struct OutputParameters {
public:
  explicit OutputParameters(const std::string _path = DSC_CONFIG_GET("global.datadir", "data"));
  virtual ~OutputParameters() {}
private:
  std::string my_prefix_;
  const std::string my_path_;
public:
  void set_prefix(std::string my_prefix);
  //! path where the data is stored
  std::string path() const;
private:
  //! base of file name for data file
  std::string prefix() const;

  /** format of output:
   *  0; // GRAPE (lossless format)
   *  1;       // VTK
   *  2; // VTK vertex data
   *  3; // gnuplot
  **/
  int outputformat() const;

  //! return true if all data should be written to a spearate path per rank
  virtual bool separateRankPath() const;

  virtual std::string macroGridName(const int dim) const;
public:
  std::string fullpath(const std::string function);

private:
  //! to avoid confusion path is only changeable in ctor
  void set_path(std::string my_path) = delete;
};

} // namespace Multiscale
} // namespace Dune

#endif // OUTPUTPARAMETER_HH
