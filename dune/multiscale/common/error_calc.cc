#include "error_calc.hh"

#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/l2error.hh>

#include <dune/multiscale/problems/elliptic/selector.hh>
#include <iostream>

Dune::Multiscale::ErrorCalculator::ErrorCalculator(const CommonTraits::DiscreteFunctionType *msfem_solution,
                                                   const CommonTraits::DiscreteFunctionType *fem_solution)
    : msfem_solution_(msfem_solution)
    , fem_solution_(fem_solution)
{}

void Dune::Multiscale::ErrorCalculator::print(std::ostream &out)
{
    assert(msfem_solution_ || fem_solution_);
    out << std::endl << "The L2 errors:" << std::endl << std::endl;

    auto gridPart = fem_solution_ ? fem_solution_->gridPart()
                                  : msfem_solution_->gridPart();
    Dune::H1Norm< CommonTraits::GridPartType > h1norm(gridPart);
    Dune::L2Error< typename CommonTraits::DiscreteFunctionType > l2error;

    //! ----------------- compute L2- and H1- errors -------------------
    if (Problem::ModelProblemData::has_exact_solution)
    {

      const CommonTraits::ExactSolutionType u;
      const int experimentally_determined_maximum_order_for_GridFunctionAdapter_bullshit = 6;
      const Dune::GridFunctionAdapter<CommonTraits::ExactSolutionType, CommonTraits::GridPartType>
          u_disc("", u, gridPart, experimentally_determined_maximum_order_for_GridFunctionAdapter_bullshit);

      if (msfem_solution_)
      {
          CommonTraits::RangeType msfem_error = l2error.norm(u, *msfem_solution_ );
          out << "|| u_msfem - u_exact ||_L2 =  " << msfem_error << std::endl;

          CommonTraits::RangeType h1_msfem_error = h1norm.distance(u_disc, *msfem_solution_);
          out << "|| u_msfem - u_exact ||_H1 =  " << h1_msfem_error << std::endl << std::endl;
      }

      if (fem_solution_)
      {
        CommonTraits::RangeType fem_error = l2error.norm(u, *fem_solution_);
        out << "|| u_fem - u_exact ||_L2 =  " << fem_error << std::endl;

        CommonTraits::RangeType h1_fem_error = h1norm.distance(u_disc, *fem_solution_);
        out << "|| u_fem - u_exact ||_H1 =  " << h1_fem_error << std::endl << std::endl;
      }
    }
    if ( msfem_solution_ && fem_solution_) {
      CommonTraits::RangeType approx_msfem_error = l2error.norm2< 2* CommonTraits::DiscreteFunctionSpaceType::polynomialOrder + 2 >(*fem_solution_,
                                                                                                        *msfem_solution_);
      out << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std::endl;

      CommonTraits::RangeType h1_approx_msfem_error = h1norm.distance(*fem_solution_, *msfem_solution_);
      out << "|| u_msfem - u_fem ||_H1 =  " << h1_approx_msfem_error << std::endl << std::endl;
    }
}
