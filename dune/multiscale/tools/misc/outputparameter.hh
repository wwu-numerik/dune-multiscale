// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef OUTPUTPARAMETER_HH
#define OUTPUTPARAMETER_HH

#include <string>

#include <dune/fem/io/file/dataoutput.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/filesystem.hh>

namespace Dune {
//! define output parameters for \ref Dune::DataOutput
struct myDataOutputParameters
  : public Dune::DataOutputParameters
{
public:
  explicit myDataOutputParameters(const std::string _path
                                  = DSC_CONFIG_GET("global.datadir", "data"))
    : my_prefix_("solutions")
    , my_path_(_path)
  {
    Dune::Stuff::Common::testCreateDirectory(my_path_);
  }

  std::string my_prefix_;
  const std::string my_path_;

  void set_prefix(std::string my_prefix) {
    my_prefix_ = my_prefix;
    Dune::Stuff::Common::testCreateDirectory(my_prefix_);
    // std :: cout << "Set prefix. my_prefix_ = " << my_prefix_ << std :: endl;
  }

  //! base of file name for data file
  std::string prefix() const {
      return my_prefix_;
  }

  //! path where the data is stored
  std::string path() const {
      return my_path_;
  }

  /** format of output:
   *  0; // GRAPE (lossless format)
   *  1;       // VTK
   *  2; // VTK vertex data
   *  3; // gnuplot
  **/
  int outputformat() const {
    // return 0; // GRAPE (lossless format)
    return 1;       // VTK
    // return 2; // VTK vertex data
    // return 3; // gnuplot
  }

private:
  //! to avoid confusion path is only changeable in ctor
  void set_path(std::string my_path) = delete;
};

} //namespace Dune
#endif // OUTPUTPARAMETER_HH
