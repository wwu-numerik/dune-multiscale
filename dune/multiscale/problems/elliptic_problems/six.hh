#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SIX
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SIX

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

// ! USE THIS ONE FOR MSFEM TESTS! PURELY LINEAR ELLIPTIC!

// !############################## Elliptic Problem 6 ###################################


// eps = 0.001 => H Ref = 16
// eps = 0.002 => H Ref = 14
// eps = 0.004 => H Ref = 12
// eps = 0.008 => H Ref = 10
// eps = 0.016 => H Ref = 8
// eps = 0.032 => H Ref = 6
// eps = 0.064 => H Ref = 4
// eps = 0.128 => H Ref = 2
// eps = 0.256 => H Ref = 0

// NOTE that (delta/epsilon_est) needs to be a positive integer!

// is an exact solution available?
// #define EXACTSOLUTION_AVAILABLE

// Note that in the following, 'Imp' abbreviates 'Implementation'
namespace Problem {
namespace Six {
// description see below 0.05
CONSTANTSFUNCTION(0.01, 0.01, 0.01)

// model problem information
struct ModelProblemData
  : public IModelProblemData
{
  ModelProblemData(const std::string filename = "no_name")
    : IModelProblemData(constants(), filename) {
  }

  inline int get_Number_of_Model_Problem() const {
    return 6;
  }

  //! \copydoc IModelProblemData::getMacroGridFile()
  inline std::string getMacroGridFile() const {
    return("../dune/multiscale/grids/macro_grids/elliptic/cube_one.dgf");
  }

  //! \copydoc IModelProblemData::getRefinementLevelReferenceProblem()
  inline int getRefinementLevelReferenceProblem() const {
    return 10;
  }
};

CONSTANTFUNCTION(FirstSource, 1.0)
NULLFUNCTION(SecondSource)

template< class FunctionSpaceImp >
class Diffusion
  : public Dune::Fem::Function< FunctionSpaceImp, Diffusion< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef Diffusion< FunctionSpaceType >                     ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
    Diffusion(){}

  // in the linear setting, use the structure
  // A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
  // the usage of an evaluate method with "evaluate ( i, j, x, y, z)" should be avoided
  // use "evaluate ( i, x, y, z)" instead and return RangeType-vector.

  // instantiate all possible cases of the evaluate-method:

  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& x,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    const auto coeff = constants().coefficients(x);
    flux[0][0] = coeff.first * gradient[0][0];
    flux[0][1] = coeff.second * gradient[0][1];
  } // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const {
    const auto coeff = constants().coefficients(x);
    flux[0][0] = coeff.first * direction_gradient[0][0];
    flux[0][1] = coeff.second * direction_gradient[0][1];
  } // jacobianDiffusiveFlux

// deprecated
// ! defaults (not to be used):
// z_i = A^{\epsilon}_i(x,vec)
// instantiate all possible cases of the evaluate-method:
  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       RangeType& z) const {
    if (i == j)
    { z = 1.01 + cos(2.0 * M_PI * x[0]); }       // !x[ 0 ] + x[ 1 ];} //!cos( 2.0 * M_PI * x[0] ); }
    else
    { z = 0.0; }
  }

  inline void evaluate(const int /*i*/,
                       const int /*j*/,
                       const DomainType& /*x*/,
                       const DomainType& /*y*/,
                       RangeType& /*z*/) const {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  } // evaluate

  inline void evaluate(const int /*i*/,
                       const int /*j*/,
                       const DomainType& /*x*/,
                       const TimeType& /*time*/,
                       RangeType& /*z*/) const {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  } // evaluate

  // dummy implementation
  inline void evaluate(const DomainType& /*x*/,
                       RangeType& /*y*/) const {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  } // evaluate
};

template< class FunctionSpaceImp, class FieldMatrixImp >
class HomDiffusion
  : public Dune::Fem::Function< FunctionSpaceImp, HomDiffusion< FunctionSpaceImp, FieldMatrixImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;
  typedef FieldMatrixImp   FieldMatrixType;

private:
  typedef HomDiffusion< FunctionSpaceType, FieldMatrixType > ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  const FieldMatrixType& A_hom_;

public:
  inline explicit HomDiffusion(const FieldMatrixType& A_hom)
    : A_hom_(A_hom)
  {}

  // in the linear setting, use the structure
  // A^{\epsilon}_i(x,\xi) = A^{\epsilon}_{i1}(x) \xi_1 + A^{\epsilon}_{i2}(x) \xi_2
  // the usage of an evaluate method with "evaluate ( i, j, x, y, z)" should be avoided
  // use "evaluate ( i, x, y, z)" instead and return RangeType-vector.

  // instantiate all possible cases of the evaluate-method:

  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& /*x*/,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    flux[0][0] = A_hom_[0][0] * gradient[0][0] + A_hom_[0][1] * gradient[0][1];
    flux[0][1] = A_hom_[1][0] * gradient[0][0] + A_hom_[1][1] * gradient[0][1];
  }

  /** the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
   * "nabla w", i.e.
   * jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:
   * jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient 
  **/
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

CONSTANTFUNCTION(MassTerm,  0.00001)
NULLFUNCTION(DefaultDummyFunction)
NULLFUNCTION(ExactSolution)

} // namespace Six {
}


#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SIX
