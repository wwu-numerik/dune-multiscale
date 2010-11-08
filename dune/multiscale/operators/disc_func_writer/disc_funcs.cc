/**
 *  \file   disc_funcs.cc
 *
 *  \brief  brief
 **/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/fem/misc/mpimanager.hh> // An initializer of MPI
#include <dune/common/exceptions.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh> // for the grid
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/dgspace/dgadaptiveleafgridpart.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction

#include <iostream>

#include "discretefunctionwriter.hh"

template < class Function >
void addScalarToFunc( Function& f, double sc )
{
    typedef typename Function::DofIteratorType
        DofIteratorType;
    DofIteratorType it = f.dbegin();
    for ( ; it != f.dend(); ++it )
        *it += sc;
    return;
}

template < class Stream, class DiscFunc >
void oneLinePrint( Stream& stream, const DiscFunc& func )
{
    typedef typename DiscFunc::ConstDofIteratorType
        DofIteratorType;
    DofIteratorType it = func.dbegin();
    stream << "\n" << func.name() << ": [ ";
    for ( ; it != func.dend(); ++it )
        stream << std::setw(5) << *it << "  ";

    stream << " ] " << std::endl;
}

/**
 *  \brief  main function
 *
 *  \attention  attention
 *
 *  \param  argc
 *          number of arguments from command line
 *  \param  argv
 *          array of arguments from command line
 **/
int main( int argc, char** argv )
{
    try
    {
        Dune::MPIManager::initialize(argc, argv);

        Dune::GridPtr< GridType > gridPtr( "grid_2d.dgf" );
        typedef Dune::LeafGridPart< GridType >
            GridPartType;
        GridPartType gridPart( *gridPtr );
        const int gridDim = GridType::dimensionworld;
        const int polOrder = POLORDER;

        typedef Dune::FunctionSpace< double, double, gridDim, gridDim >
            FunctionSpaceType;
        typedef Dune::DiscontinuousGalerkinSpace<   FunctionSpaceType,
                                                    GridPartType,
                                                    polOrder >
            DGSpaceType;
        typedef Dune::AdaptiveDiscreteFunction< DGSpaceType >
            DFType;
        DGSpaceType dg_space( gridPart );
        DFType func1 ( "f1", dg_space );
        DFType func2 ( "f2", dg_space );
        DFType func3 ( "f3", dg_space );
        addScalarToFunc( func1, 1 );
        addScalarToFunc( func2, 2 );
        addScalarToFunc( func3, 3 );
        bool ok = false;

        {
            DiscreteFunctionWriter dfw( "dump" );
            ok = dfw.open();

            if ( ok )
            {
                dfw.append( func1 );
                dfw.append( func2 );
                dfw.append( func3 );
                func1.clear();
                func2.clear();
                func3.clear();
            }
        }

        {
            DiscreteFunctionReader dfr( "dump" );
            ok = dfr.open();

            if ( ok )
            {
                dfr.read( 0, func1 );
                dfr.read( 1, func2 );
                dfr.read( 2, func3 );
                oneLinePrint( std::cout , func1 );
                oneLinePrint( std::cout , func2 );
                oneLinePrint( std::cout , func3 );
                func2 += func1;
                func3 -= func2;
                oneLinePrint( std::cout , func3 );
            }
        }
        return ok ? 0 : -1;
    }
    catch (Dune::Exception &e)
    {
        std::cerr << "Dune reported error: " << e << std::endl;
        return -1;
    }
    catch (...)
    {
        std::cerr << "Unknown exception thrown!" << std::endl;
        return -1;
    }
}

