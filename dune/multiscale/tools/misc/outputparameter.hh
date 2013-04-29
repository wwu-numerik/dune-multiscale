// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef OUTPUTPARAMETER_HH
#define OUTPUTPARAMETER_HH

#include <string>

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>

namespace Dune {
namespace Multiscale {
//! define output parameters for \ref Dune::DataOutput
struct OutputParameters
  : public Dune::DataOutputParameters
{
public:
  explicit OutputParameters(const std::string _path
                                  = DSC_CONFIG_GET("global.datadir", "data"));

  std::string my_prefix_;
  const std::string my_path_;

  void set_prefix(std::string my_prefix);
  //! base of file name for data file
  std::string prefix() const;
  //! path where the data is stored
  std::string path() const;

  /** format of output:
   *  0; // GRAPE (lossless format)
   *  1;       // VTK
   *  2; // VTK vertex data
   *  3; // gnuplot
  **/
  int outputformat() const;

private:
  //! to avoid confusion path is only changeable in ctor
  void set_path(std::string my_path) = delete;
};

} //namespace Multiscale
} //namespace Dune

#endif // OUTPUTPARAMETER_HH
