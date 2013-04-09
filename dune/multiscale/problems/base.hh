// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS_PROBLEMS_BASE_HH
#define DUNE_MS_PROBLEMS_BASE_HH

#include <string>
#include <dune/multiscale/problems/constants.hh>
#include <dune/stuff/fem/functions/analytical.hh>

namespace Problem {
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
 *   evaluate  \grad u^{\epsilon}( x )  --> evaluateJacobian


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
class IModelProblemData
{
protected:
  const Constants constants_;

public:

  static const bool has_exact_solution = false;
  //! Constructor for ModelProblemData
  inline IModelProblemData(const Constants constants)
    : constants_(constants)
  {}
  virtual ~IModelProblemData()
  {}

  /**
   * @brief getMacroGridFile returns a path to a Dune::Grid loadable file (dgf)
   * @return macroGridName is set to said path
   * \todo paths need to be relative to binary
   */
  virtual std::string getMacroGridFile() const = 0;

  // are the coefficients periodic? (e.g. A=A(x/eps))
  // this method is only relevant if you want to use a standard homogenizer
  virtual bool problemIsPeriodic() const = 0;

  // does the problem allow a stochastic perturbation of the coefficients?
  virtual bool problemAllowsStochastics() const = 0;

};

// linear constant diffusion operator that can be filled with given values
// (e.g. with the pre-computed values of a homogenzid matrix)
template< class FunctionSpaceImp, class FieldMatrixImp >
class ConstantDiffusionMatrix
  : public Dune::Fem::Function< FunctionSpaceImp, ConstantDiffusionMatrix< FunctionSpaceImp, FieldMatrixImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;
  typedef FieldMatrixImp   FieldMatrixType;

private:
  typedef ConstantDiffusionMatrix< FunctionSpaceType, FieldMatrixType > ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  const FieldMatrixType& A_values_;

public:
  inline explicit ConstantDiffusionMatrix(const FieldMatrixType& A_given_values)
    : A_values_(A_given_values)
  {}

  // (diffusive) flux = A( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& /*x*/,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    flux[0][0] = A_values_[0][0] * gradient[0][0] + A_values_[0][1] * gradient[0][1];
    flux[0][1] = A_values_[1][0] * gradient[0][0] + A_values_[1][1] * gradient[0][1];
  }

  // should not be required, since we are in a fully linear setting
  void jacobianDiffusiveFlux(const DomainType&,
                             const JacobianRangeType&,
                             const JacobianRangeType&,
                             JacobianRangeType&) const {
    DUNE_THROW(Dune::NotImplemented,"");
  } // jacobianDiffusiveFlux

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  }
};

} //! @} namespace Problem

#endif // DUNE_MS_PROBLEMS_BASE_HH
