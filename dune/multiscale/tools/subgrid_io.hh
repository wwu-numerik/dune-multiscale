#ifndef DUNE_MULTISCALE_SUBGRID_IO_HH
#define DUNE_MULTISCALE_SUBGRID_IO_HH

/** this file contains specializations of write/read functions targetd for dune-subgrid
 **/

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/utility/grapedataioformattypes.hh>
#include <dune/grid/yaspgrid.hh>
//#ifdef HAVE_ALUGRID
  #include <dune/grid/alugrid.hh>
//#endif
#include <dune/grid/sgrid.hh>
#include <dune/subgrid/subgrid.hh>
#include <boost/filesystem/fstream.hpp>

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>

namespace Dune {

    template <typename HostGridType>
    bool writeHostGrid(HostGridType& hostgrid, std::string filename);
    template <typename HostGridType>
    bool readHostGrid(HostGridType& hostgrid, std::string filename);

    template <class HostgridType, class SubgridType>
    std::string subgridKeygen(const HostgridType& hostgrid, const SubgridType& subgrid, const int subgrid_idx)
    {
        boost::format key("%s_%s_s%d");
        key % hostgrid.name() % subgrid.name() % subgrid_idx;
        return key.str();
    }

    template<class AlugridType>
    bool writeAlugrid(const AlugridType& grid, const std::string filename)
    {
        if (boost::filesystem::exists(filename))
            return true;
        grid.template writeGrid<GrapeIOFileFormatType::xdr>(filename, 0.0f);
        return true;
    }

    template<class AlugridType>
    bool readAlugrid(AlugridType& grid, const std::string filename)
    {
        if (grid.size(0))
            return true;
        double dummy;
        grid.template readGrid<GrapeIOFileFormatType::xdr>(filename, dummy);
        return true;
    }

#define ALUGRID_IO_FUNCTION_PAIR(classname,dim) \
    template<> bool writeHostGrid(classname<dim,dim>& hostgrid, std::string filename) \
    { return writeAlugrid(hostgrid, filename); }\
    \
    template<> bool readHostGrid(classname<dim,dim>& hostgrid, std::string filename) \
    { return readAlugrid(hostgrid, filename); }

    ALUGRID_IO_FUNCTION_PAIR(ALUSimplexGrid,2)
    ALUGRID_IO_FUNCTION_PAIR(ALUConformGrid,2)
    ALUGRID_IO_FUNCTION_PAIR(ALUCubeGrid,2)
    ALUGRID_IO_FUNCTION_PAIR(ALUSimplexGrid,3)
    ALUGRID_IO_FUNCTION_PAIR(ALUCubeGrid,3)
#undef ALUGRID_IO_FUNCTION_PAIR

//#define ALUGRID_IO_FUNCTION_PAIR(classname) \
//    template<class E, class C> bool writeHostGrid(classname<E,C>& hostgrid, std::string filename) \
//    { return writeAlugrid(hostgrid, filename); }\
//    \
//    template<class E, class C> bool readHostGrid(classname<E,C>& hostgrid, std::string filename) \
//    { return readAlugrid(hostgrid, filename); }

//    ALUGRID_IO_FUNCTION_PAIR(ALU3dGrid)
//    ALUGRID_IO_FUNCTION_PAIR(ALU2dGrid)
//#undef ALUGRID_IO_FUNCTION_PAIR

}

#endif // SUBGRID_IO_HH
