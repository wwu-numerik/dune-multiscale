// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS_PROBLEMS_BASE_HH
#define DUNE_MS_PROBLEMS_BASE_HH


#include <dune/multiscale/common/traits.hh>
#include <dune/multiscale/msfem/msfem_traits.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/fem/functions/analytical.hh>
#include <dune/stuff/functions/global.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/global.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <memory>
#include <string>

#include "config.h"

namespace Dune {
namespace Multiscale {
namespace Problem {

struct DiffusionBase {

  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  virtual ~DiffusionBase() {}

  //! in the linear setting, use the structure
  //! A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
  //! (diffusive) flux = A^{\epsilon}( x , direction )
  //! (typically direction is some 'gradient_of_a_function')
  virtual void diffusiveFlux(const DomainType& x, const JacobianRangeType& direction,
                             JacobianRangeType& flux) const = 0;

  //! the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  //! "nabla w", i.e.
  //! jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:
  //! jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  virtual void jacobianDiffusiveFlux(const DomainType& x, const JacobianRangeType& /*position_gradient*/,
                                     const JacobianRangeType& direction_gradient, JacobianRangeType& flux) const = 0;
};

struct LowerOrderBase : public Dune::Multiscale::CommonTraits::FunctionBaseType {
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef double TimeType;


  virtual void evaluate(const DomainType& x, RangeType& ret) const { evaluate(x, TimeType(0), ret); }
  virtual void evaluate(const DomainType& x, const TimeType& time, RangeType& y) const = 0;
  virtual void evaluate(const DomainType& x, const RangeType& position, const JacobianRangeType& direction_gradient,
                        RangeType& y) const = 0;
  virtual void position_derivative(const DomainType& x, const RangeType& position,
                                   const JacobianRangeType& direction_gradient, RangeType& y) const = 0;
  virtual void direction_derivative(const DomainType& x, const RangeType& position,
                                    const JacobianRangeType& direction_gradient, JacobianRangeType& y) const = 0;
};

struct ZeroLowerOrder : public LowerOrderBase {
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef double TimeType;

  virtual void evaluate(const DomainType&, RangeType& y) const DS_FINAL { y = RangeType(0); }
  virtual void evaluate(const DomainType&, const TimeType&, RangeType& y) const DS_FINAL{ y = RangeType(0); }
  virtual void evaluate(const DomainType&, const RangeType&, const JacobianRangeType&, RangeType& y) const DS_FINAL {
    y = RangeType(0);
  }
  virtual void position_derivative(const DomainType&, const RangeType&, const JacobianRangeType&, RangeType& y) const DS_FINAL {
    y = RangeType(0);
  }
  virtual void direction_derivative(const DomainType&, const RangeType&, const JacobianRangeType&,
                                    JacobianRangeType& y) const DS_FINAL {
    y = JacobianRangeType(0);
  }
};

class DirichletDataBase : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef DomainFieldType TimeType;

  virtual void evaluate(const DomainType& x, RangeType& y) const = 0;
  virtual void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const = 0;

};

class ZeroDirichletData : public DirichletDataBase {
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef DomainFieldType TimeType;

  virtual void evaluate(const DomainType& /*x*/, RangeType& y) const DS_FINAL { y = RangeType(0.0); }
  virtual void evaluate(const DomainType& /*x*/, const TimeType& /*time*/, RangeType& y) const DS_FINAL { y = RangeType(0.0); }
};

class NeumannDataBase : public Dune::Multiscale::CommonTraits::FunctionBaseType {
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef DomainFieldType TimeType;

  virtual void evaluate(const DomainType& x, RangeType& y) const = 0;
  virtual void evaluate(const DomainType& x, const TimeType& /*time*/, RangeType& y) const = 0;

};

class ZeroNeumannData : public NeumannDataBase {
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef DomainFieldType TimeType;

  virtual void evaluate(const DomainType& /*x*/, RangeType& y) const DS_FINAL { y = RangeType(0.0); }
  virtual void evaluate(const DomainType& /*x*/, const TimeType& /*time*/, RangeType& y) const DS_FINAL { y = RangeType(0.0); }
};

/**
 * \addtogroup Problem
 * @{
 *
 * in general we regard problems of the following type:

 * - div ( A^{\epsilon} (x,\nabla u^{\epsilon}) ) + m^{\epsilon} u^{\epsilon} = f - div G

 * Here we have:
 * u^{\epsilon} = exact solution of the problem
 * A^{\epsilon} = diffusion matrix (or tensor), e.g. with the structure A^{\epsilon}(x) = A(x,\frac{x}{\epsilon})
 * m^{\epsilon} = a mass term (or reaction term), e.g. with the structure m^{\epsilon}(x) = m(x,\frac{x}{\epsilon})
 * f = first source term with the structure f = f(x) (=> no micro-scale dependency)
 * G = second source term with the structure G = G(x) (=> no micro-scale dependency).
 * (Note that 'G' is directly implemented! We do not implement '- div G'!)

 * A^{\epsilon} is can be a monotone operator (=> use HMM, the MsFEM is not implemented for nonlinear problems)


 * ! we use the following class names:


 * class ExactSolution -> describes u^{\epsilon}
 * methods:
 *   evaluate  u^{\epsilon}( x )        --> evaluate
 *   evaluate  \grad u^{\epsilon}( x )  --> jacobian


 * class Diffusion -> describes A^{\epsilon}
 * methods:
 *   evaluate A^{\epsilon}( x , direction )            --> diffusiveFlux
 *   evaluate DA^{\epsilon}( x , position ) direction  --> jacobianDiffusiveFlux


 * class MassTerm -> describes m^{\epsilon}
 * methods:
 *   evaluate m^{\epsilon}( x )         --> evaluate
 *

 * class FirstSource -> describes f
 * methods:
 *   evaluate f( x )                    --> evaluate


 * class SecondSource -> describes G
 * methods:
 *   evaluate G( x )                    --> evaluate


 ! See 'elliptic_problems/example.hh' for details


 * The mass (or reaction) term m^{\epsilon} is given by:
 * ! m^{\epsilon} := \epsilon
 * Since \epsilon tends to zero, we may say that we do not have a real mass term for our problem. It is a simple
 * condition to fix the solution which is only unique up to a constant. In fact we still approximate the solution of the
 * problem without mass.

 * The first source term f is given by:
 * ! f(x) := ****
 * since the second source is zero, f will form the right hand side (RHS) of our discrete problem

 * The second source term G is constantly zero:
 * ! G(x) := 0

 * !FirstSource defines the right hand side (RHS) of the governing problem (i.e. it defines 'f').
 * The value of the right hand side (i.e. the value of 'f') at 'x' is accessed by the method 'evaluate'. That means 'y
 * := f(x)' and 'y' is returned. It is only important that 'RHSFunction' knows the function space ('FuncSpace') that it
 * is part from. (f \in FunctionSpace)
 *
 *
**/
class IModelProblemData {
protected:
  const Constants constants_;
  typedef CommonTraits::GridType::LeafGridView View;
  typedef Dune::Stuff::GridboundaryInterface<typename View::Intersection> BoundaryInfoType;
  typedef MsFEM::MsFEMTraits::LocalGridType::LeafGridView SubView;
  typedef Dune::Stuff::GridboundaryInterface<typename SubView::Intersection> SubBoundaryInfoType;

public:

  //! Constructor for ModelProblemData
  inline IModelProblemData(const Constants constants) : constants_(constants) {}
  virtual ~IModelProblemData() {}

  /**
   * @brief getMacroGridFile returns a path to a Dune::Grid loadable file (dgf)
   * @return macroGridName is set to said path
   * \todo paths need to be relative to binary
   */
  virtual std::string getMacroGridFile() const = 0;

  //! are the coefficients periodic? (e.g. A=A(x/eps))
  //! this method is only relevant if you want to use a standard homogenizer
  virtual bool problemIsPeriodic() const = 0;

  //! does the problem allow a stochastic perturbation of the coefficients?
  virtual bool problemAllowsStochastics() const = 0;

  //! does the problem implement an exact solution?
  virtual bool hasExactSolution() const { return false; }

  //! is the diffusion matrix symmetric?
  virtual bool symmetricDiffusion() const { return true; }

  //! linear/nonlinear toggle
  virtual bool linear() const { return true; }

  virtual std::unique_ptr<BoundaryInfoType> boundaryInfo() const {
    return DSC::make_unique<DS::GridboundaryAllDirichlet<typename View::Intersection>>();
  }
  virtual std::unique_ptr<SubBoundaryInfoType> subBoundaryInfo() const {
    return DSC::make_unique<DS::GridboundaryAllDirichlet<typename SubView::Intersection>>();
  }
};

} //! @} namespace Problem
} // namespace Multiscale {
} // namespace Dune {

namespace DMP = Dune::Multiscale::Problem;

#define MSCONSTANTFUNCTION(classname, constant)                                                                        \
  class classname : public Dune::Multiscale::CommonTraits::ConstantFunctionBaseType {                                  \
  public:                                                                                                              \
    classname() : Dune::Multiscale::CommonTraits::ConstantFunctionBaseType(constant) {}                                \
  };

#define MSNULLFUNCTION(classname)                                                                                      \
  class classname : public Dune::Multiscale::CommonTraits::ConstantFunctionBaseType {                                  \
  public:                                                                                                              \
    classname() : Dune::Multiscale::CommonTraits::ConstantFunctionBaseType(0.0) {}                                     \
    classname(const Dune::Multiscale::CommonTraits::FunctionSpaceType&)                                                \
      : Dune::Multiscale::CommonTraits::ConstantFunctionBaseType(0.0) {}                                               \
  };

#endif // DUNE_MS_PROBLEMS_BASE_HH
