#include "selector.hh"

//for i in $(ls *hh) ; do echo \#include \"${i}\" ; done
#include "elliptic/eight.hh"
#include "elliptic/eleven.hh"
#include "elliptic/example.hh"
#include "elliptic/five.hh"
#include "elliptic/four.hh"
#include "elliptic/nine.hh"
#include "elliptic/one.hh"
#include "elliptic/seven.hh"
#include "elliptic/six.hh"
#include "elliptic/ten.hh"
#include "elliptic/three.hh"
#include "elliptic/toy.hh"
#include "elliptic/two.hh"

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


std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> Dune::Multiscale::Problem::getMassTerm()
{
    return DSC::make_unique<Dune::Multiscale::Problem::PROBLEM_NAME::MassTerm>();
}


std::unique_ptr<Dune::Multiscale::CommonTraits::FunctionBaseType> Dune::Multiscale::Problem::getDefaultDummyFunction()
{
    return DSC::make_unique<Dune::Multiscale::Problem::PROBLEM_NAME::DefaultDummyFunction>();
}


std::unique_ptr<Dune::Multiscale::CommonTraits::ModelProblemDataType> Dune::Multiscale::Problem::getModelData()
{
    return DSC::make_unique<Dune::Multiscale::Problem::PROBLEM_NAME::ModelProblemData>();
}


std::unique_ptr<const Dune::Multiscale::CommonTraits::LowerOrderTermType> Dune::Multiscale::Problem::getLowerOrderTerm()
{
    return DSC::make_unique<const Dune::Multiscale::Problem::PROBLEM_NAME::LowerOrderTerm>();
}


std::unique_ptr<Dune::Multiscale::CommonTraits::DiffusionType> Dune::Multiscale::Problem::getDiffusion()
{
    return DSC::make_unique<Dune::Multiscale::Problem::PROBLEM_NAME::Diffusion>();
}


std::string Dune::Multiscale::Problem::name()
{
  return "Dummy Dune::Multiscale::Problem::name()";
}
