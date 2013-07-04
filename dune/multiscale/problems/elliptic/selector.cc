#include "selector.hh"


std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> Dune::Multiscale::Problem::getFirstSource() {
    return DSC::make_unique<Dune::Multiscale::Problem::PROBLEM_NAME::FirstSource>();
}


std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> Dune::Multiscale::Problem::getSecondSource()
{
    return DSC::make_unique<Dune::Multiscale::Problem::PROBLEM_NAME::SecondSource>();
}


std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> Dune::Multiscale::Problem::getExactSolution()
{
    return DSC::make_unique<Dune::Multiscale::Problem::PROBLEM_NAME::ExactSolution>();
}
