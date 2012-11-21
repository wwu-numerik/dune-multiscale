#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FOUR
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FOUR

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

//!############################## Elliptic Problem 4 ###################################


// Note that in the following, 'Imp' abbreviates 'Implementation'
namespace Problem {
namespace Four {
// description see below 0.05
CONSTANTSFUNCTION( 0.05 )
// NOTE that (delta/epsilon_est) needs to be a positive integer!

// model problem information
struct ModelProblemData
  : public IModelProblemData
{
  ModelProblemData()
    : IModelProblemData(constants()) {
  }

  //! \copydoc IModelProblemData::getMacroGridFile()
  inline std::string getMacroGridFile() const {
    return("../dune/multiscale/grids/macro_grids/elliptic/corner_singularity.dgf");
  }

  //! \copydoc IModelProblemData::getRefinementLevelReferenceProblem()
  inline int getRefinementLevelReferenceProblem() const {
    return 20;
  }
};

template< class FunctionSpaceImp >
class FirstSource
  : public Dune::Fem::Function< FunctionSpaceImp, FirstSource< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef FirstSource< FunctionSpaceType >                   ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  static const int dimDomain = DomainType::dimension;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  FirstSource(){}

  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    // Krei um einspringende Ecke (Kreis um (0.5,0.5) )
    double distance = sqrt( pow(x[0] - 0.5, 2.0) + pow(x[1] - 0.5, 2.0) );

    if (distance < 0.2)
    { y = 1.0; } else
    { y = 0.1; }
  } // evaluate

  inline void evaluate(const DomainType& x,
                       const TimeType& /*time*/,
                       RangeType& y) const {
    evaluate(x, y);
  }
};

  /** \brief default class for the second source term G.
   * Realization: set G(x) = 0: **/
  NULLFUNCTION(SecondSource)

// the (non-linear) diffusion operator A^{\epsilon}(x,\xi)
// A^{\epsilon} : R^d -> R^d

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
    if (constants().get("linear", true)) {
      flux[0][0] = coeff.first * gradient[0][0];
      flux[0][1] = coeff.second * gradient[0][1];
    } else {
      flux[0][0] = coeff.first  * ( gradient[0][0] + ( (1.0 / 3.0) * pow(gradient[0][0], 3.0) ) );
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

    if (constants().get("linear", true)) {
      flux[0][0] = coeff.first * direction_gradient[0][0];
      flux[0][1] = coeff.second * direction_gradient[0][1];
    } else {
      flux[0][0] = coeff.first * direction_gradient[0][0]
                   * ( 1.0 + pow(position_gradient[0][0], 2.0) );
      flux[0][1] = coeff.second * direction_gradient[0][1]
                   * ( 1.0 + pow(position_gradient[0][1], 2.0) );
    }
  } // jacobianDiffusiveFlux

  /** \deprecated throws Dune::NotImplemented exception **/
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
    DUNE_THROW(Dune::NotImplemented,"No homogenization available");
    if (constants().get("linear", true)) {
      flux[0][0] = A_hom_[0][0] * gradient[0][0] + A_hom_[0][1] * gradient[0][1];
      flux[0][1] = A_hom_[1][0] * gradient[0][0] + A_hom_[1][1] * gradient[0][1];
    } else {
      flux[0][0] = A_hom_[0][0] * gradient[0][0] + A_hom_[0][1] * gradient[0][1];
      flux[0][1] = A_hom_[1][0] * gradient[0][0] + A_hom_[1][1] * gradient[0][1];
      //! TODO one of the the above is in the wrong branch
      DUNE_THROW(Dune::NotImplemented,"Nonlinear example not yet implemented.");
    }
  } // diffusiveFlux

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
} // namespace Four {
}


#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_FOUR
