// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_SUBGRID_IO_HH
#define DUNE_MULTISCALE_SUBGRID_IO_HH

/** this file contains specializations of write/read functions targetd for dune-subgrid
 **/

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#else
  #include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/utility/grapedataioformattypes.hh>
#include <dune/grid/yaspgrid.hh>
#ifdef HAVE_ALUGRID
  #include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/spgrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/subgrid/subgrid.hh>
#include <boost/filesystem/fstream.hpp>

#include <dune/stuff/aliases.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/filesystem.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>

namespace Dune {

    template <typename HostGridType>
    static bool writeHostGrid(HostGridType& hostgrid, std::string filename);
    template <typename HostGridType>
    static bool readHostGrid(HostGridType& hostgrid, std::string filename);

    template <class HostgridType, class SubgridType>
    std::string subgridKeygen(const HostgridType& /*hostgrid*/, const SubgridType& subgrid, const int subgrid_idx)
    {
        boost::format key("%s_%s_s%d");
        key % DSC::Typename<HostgridType>::value() % subgrid.name() % subgrid_idx;
        return key.str();
    }

    template<class HostgridType>
    bool writeHostgridCommon(const HostgridType& grid, const std::string filename)
    {
        if (boost::filesystem::exists(filename))
            return true;
        grid.template writeGrid<GrapeIOFileFormatType::xdr>(filename, 0.0f);
        return true;
    }

    template<class HostgridType>
    bool readHostgridCommon(HostgridType& grid, const std::string filename)
    {
        if (grid.size(0))
            return true;
        double dummy;
        grid.template readGrid<GrapeIOFileFormatType::xdr>(filename, dummy);
        return true;
    }


#define HOSTGRID_IO_FUNCTION_PAIR(classname,dim) \
    template<> bool writeHostGrid(classname<dim,dim>& hostgrid, std::string filename) \
    { return writeHostgridCommon(hostgrid, filename); }\
    \
    template<> bool readHostGrid(classname<dim,dim>& hostgrid, std::string filename) \
    { return readHostgridCommon(hostgrid, filename); }

    //! careful, this only works when using grid selector
    #if USED_ALBERTAGRID_GRIDTYPE
        HOSTGRID_IO_FUNCTION_PAIR(AlbertaGrid,2)
    #elif USED_SPGRID_GRIDTYPE
        template<> bool writeHostGrid(typename GridSelector::GridType& hostgrid, std::string filename)
        { return writeHostgridCommon(hostgrid, filename); }
        \
        template<> bool readHostGrid(typename GridSelector::GridType& hostgrid, std::string filename)
        { return readHostgridCommon(hostgrid, filename); }
    #elif USED_YASPGRID_GRIDTYPE
            template<> bool writeHostGrid(typename GridSelector::GridType& hostgrid, std::string filename)
        { /*YASPGrid can't be written to disk*/ return false; }
        \
        template<> bool readHostGrid(typename GridSelector::GridType& hostgrid, std::string filename)
        { /*YASPGrid can't be read from disk*/ return false; }
    #elif USED_ALUGRID_SIMPLEX_GRIDTYPE
        HOSTGRID_IO_FUNCTION_PAIR(ALUSimplexGrid,2)
        HOSTGRID_IO_FUNCTION_PAIR(ALUSimplexGrid,3)
    #elif USED_ALUGRID_CONFORM_GRIDTYPE
        HOSTGRID_IO_FUNCTION_PAIR(ALUConformGrid,2)
    #elif USED_ALUGRID_CUBE_GRIDTYPE
        HOSTGRID_IO_FUNCTION_PAIR(ALUCubeGrid,2)
        HOSTGRID_IO_FUNCTION_PAIR(ALUCubeGrid,3)
    #endif
#undef HOSTGRID_IO_FUNCTION_PAIR



} // namespace Dune

#endif // SUBGRID_IO_HH
