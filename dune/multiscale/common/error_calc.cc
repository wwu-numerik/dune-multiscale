#include "error_calc.hh"

#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/stuff/common/filesystem.hh>

#include <dune/multiscale/problems/elliptic/selector.hh>
#include <iostream>

template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols >
struct TimeFunctionAdapter : public Dune::Stuff::FunctionInterface< DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols >
{
  typedef Dune::Stuff::FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
    WrappedType;

  TimeFunctionAdapter(const WrappedType& wr)
    : wrapped_(wr)
  {}

  virtual void evaluate(const typename WrappedType::DomainType& x,
                        typename WrappedType::RangeType& ret) const
  {
    wrapped_(x, ret);
  }

  virtual void evaluate(const typename WrappedType::DomainType& x,
                        const double& /*t*/,
                        typename WrappedType::RangeType& ret) const
  {
    wrapped_(x, ret);
  }

  const WrappedType& wrapped_;
};

template< class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols >
TimeFunctionAdapter<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
timefunctionAdapted(const Dune::Stuff::FunctionInterface< DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols >& wrapped)
{
  return TimeFunctionAdapter<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>(wrapped);
}

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
    Dune::Fem::H1Norm< CommonTraits::GridPartType > h1norm(gridPart);
    Dune::L2Error< typename CommonTraits::DiscreteFunctionType > l2error;

    std::map<std::string,double> csv;

    //! ----------------- compute L2- and H1- errors -------------------
    if (Problem::ModelProblemData::has_exact_solution)
    {

      auto u_ptr = Dune::Multiscale::Problem::getExactSolution();
      const auto& u = *u_ptr;
      const int experimentally_determined_maximum_order_for_GridFunctionAdapter_bullshit = 6;
      const Dune::Fem::GridFunctionAdapter<CommonTraits::ExactSolutionType, CommonTraits::GridPartType>
          u_disc("", u, gridPart, experimentally_determined_maximum_order_for_GridFunctionAdapter_bullshit);

      if (msfem_solution_)
      {
          CommonTraits::RangeType msfem_error = l2error.norm(timefunctionAdapted(u), *msfem_solution_ );
          out << "|| u_msfem - u_exact ||_L2 =  " << msfem_error << std::endl;

          CommonTraits::RangeType h1_msfem_error = h1norm.distance(u_disc, *msfem_solution_);
          out << "|| u_msfem - u_exact ||_H1 =  " << h1_msfem_error << std::endl << std::endl;

          csv["msfem_exact_L2"] = msfem_error;
          csv["msfem_exact_H1"] = h1_msfem_error;
      }

      if (fem_solution_)
      {
        CommonTraits::RangeType fem_error = l2error.norm(timefunctionAdapted(u), *fem_solution_);
        out << "|| u_fem - u_exact ||_L2 =  " << fem_error << std::endl;

        CommonTraits::RangeType h1_fem_error = h1norm.distance(u_disc, *fem_solution_);
        out << "|| u_fem - u_exact ||_H1 =  " << h1_fem_error << std::endl << std::endl;

        csv["fem_exact_L2"] = fem_error;
        csv["fem_exact_H1"] = h1_fem_error;
      }
    }
    if ( msfem_solution_ && fem_solution_) {
      CommonTraits::RangeType approx_msfem_error = l2error.norm2< 2* CommonTraits::DiscreteFunctionSpaceType::polynomialOrder + 2 >(*fem_solution_,
                                                                                                        *msfem_solution_);
      out << "|| u_msfem - u_fem ||_L2 =  " << approx_msfem_error << std::endl;

      CommonTraits::RangeType h1_approx_msfem_error = h1norm.distance(*fem_solution_, *msfem_solution_);
      out << "|| u_msfem - u_fem ||_H1 =  " << h1_approx_msfem_error << std::endl << std::endl;

      csv["msfem_fem_L2"] = approx_msfem_error;
      csv["msfem_fem_H1"] = h1_approx_msfem_error;
    }

    std::unique_ptr<boost::filesystem::ofstream>
            csvfile(DSC::make_ofstream(DSC_CONFIG_GET("global.datadir", "data")+"/errors.csv"));
    const std::string sep(",");
    for(const auto& key_val : csv)
    {
      *csvfile << key_val.first << sep;
    }
    *csvfile << std::endl;
    for(const auto& key_val : csv)
    {
      *csvfile << key_val.second << sep;
    }
    *csvfile << std::endl;
}
