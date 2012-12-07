//! MUSS NOCH AUFGERAEUMT UND UEBERARBEITET WERDEN!!!!!!!!!

#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SIX
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SIX

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

//!############################## Elliptic Problem 6 ###################################

// Note that in the following, 'Imp' abbreviates 'Implementation'
namespace Problem {
namespace Six {
CONSTANTSFUNCTION( 0.05 )

// model problem information
struct ModelProblemData
  : public IModelProblemData
{

  static const bool has_exact_solution = false;
  
  ModelProblemData()
    : IModelProblemData(constants()) {
      assert( constants_.epsilon != 0.0);
      if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
         DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");
  }

  //! \copydoc IModelProblemData::getMacroGridFile()
  inline std::string getMacroGridFile() const {
    return("../dune/multiscale/grids/macro_grids/elliptic/cube_two.dgf");
  }

  // are the coefficients periodic? (e.g. A=A(x/eps))
  // this method is only relevant if you want to use a standard homogenizer
  inline bool problemIsPeriodic() const {
    return true; // = problem is periodic
  }
  
  // does the problem allow a stochastic perturbation of the coefficients?
  inline bool problemAllowsStochastics() const {
    return true; // = problem allows stochastic perturbations
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
    const auto coeff = constants().coefficients_variant_A(x);
    if ( constants().get("linear", true) )
    {
      flux[0][0] = coeff.first * gradient[0][0];
      flux[0][1] = coeff.second * gradient[0][1];
    } else {
      flux[0][0] = coeff.first * ( gradient[0][0] + ( (1.0 / 3.0) * pow(gradient[0][0], 3.0) ) );
      flux[0][1] = coeff.second * ( gradient[0][1] + ( (1.0 / 3.0) * pow(gradient[0][1], 3.0) ) );
    }
  } // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& position_gradient,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const {
    const auto coeff = constants().coefficients_variant_A(x);
    if ( constants().get("linear", true) )
    {
      flux[0][0] = coeff.first * direction_gradient[0][0];
      flux[0][1] = coeff.second * direction_gradient[0][1];
    } else {
      flux[0][0] = coeff.first * direction_gradient[0][0]
                   * ( 1.0 + pow(position_gradient[0][0], 2.0) );
      flux[0][1] = coeff.second * direction_gradient[0][1]
                   * ( 1.0 + pow(position_gradient[0][1], 2.0) );
    }
  } // jacobianDiffusiveFlux

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  }
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
    if ( constants().get("linear", true) )
    {
      flux[0][0] = A_hom_[0][0] * gradient[0][0] + A_hom_[0][1] * gradient[0][1];
      flux[0][1] = A_hom_[1][0] * gradient[0][0] + A_hom_[1][1] * gradient[0][1];
    } else {
      flux[0][0] = A_hom_[0][0] * gradient[0][0] + A_hom_[0][1] * gradient[0][1];
      flux[0][1] = A_hom_[1][0] * gradient[0][0] + A_hom_[1][1] * gradient[0][1];
      //! TODO one of the the above is in the wrong branch
      DUNE_THROW(Dune::NotImplemented,"Nonlinear example not yet implemented.");
    }
  } // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& /*x*/,
                             const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& /*direction_gradient*/,
                             JacobianRangeType& /*flux*/) const {
    if ( constants().get("linear", true) )
    {
      DUNE_THROW(Dune::NotImplemented,"linear example not yet implemented.");
    } else {
      DUNE_THROW(Dune::NotImplemented,"Nonlinear example not yet implemented.");
    }
  } // jacobianDiffusiveFlux

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  }
};

CONSTANTFUNCTION(MassTerm,  0.00001)

//! a dummy function class for functions, vectors and matrices
NULLFUNCTION(DefaultDummyFunction)
NULLFUNCTION(ExactSolution)

} //namespace Six {
}

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SIX
