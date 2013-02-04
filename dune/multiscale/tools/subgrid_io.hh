#ifndef DUNE_MULTISCALE_SUBGRID_IO_HH
#define DUNE_MULTISCALE_SUBGRID_IO_HH

/** this file contains specializations of write/read functions targetd for dune-subgrid
 **/

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#ifdef HAVE_ALUGRID
  #include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>
#include <boost/filesystem/fstream.hpp>

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/filesystem.hh>

namespace Dune {

    template<class AlugridType>
    bool writeAlugrid(const AlugridType& grid, const std::string filename)
    {
        if (boost::filesystem::exists(filename))
            return true;
        grid.backup(boost::filesystem::ofstream(filename));
    }

    template<class AlugridType>
    bool readAlugrid(AlugridType& grid, const std::string filename)
    {
        if (grid.size(0))
            return true;
        grid.restore(boost::filesystem::ifstream(filename));
    }

#define ALUGRID_IO_FUNCTION_PAIR(classname) \
    template <int... A, typename... T> \
    bool writeHostGrid(classname < A..., T... >& hostgrid, std::string filename) \
    { return writeAlugrid(hostgrid, filename); } \
    \
    template <int ...A, typename ...T> \
    bool readHostGrid(classname < A..., T... >& hostgrid, std::string filename) \
    { return readAlugrid(hostgrid, filename); }

    ALUGRID_IO_FUNCTION_PAIR(ALUSimplexGrid)
    ALUGRID_IO_FUNCTION_PAIR(ALUConformGrid)
    ALUGRID_IO_FUNCTION_PAIR(ALUCubeGrid)
    ALUGRID_IO_FUNCTION_PAIR(ALU3dGrid)
    ALUGRID_IO_FUNCTION_PAIR(ALU2dGrid)

#undef ALUGRID_IO_FUNCTION_PAIR
}

#endif // SUBGRID_IO_HH
