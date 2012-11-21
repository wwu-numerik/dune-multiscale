#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EIGHT
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EIGHT

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

//!############################## Elliptic Problem 1 ###################################


// Note that in the following, 'Imp' abbreviates 'Implementation'
namespace Problem {
namespace Eight {
// description see below
// vorher: 0.13
CONSTANTSFUNCTION( 0.0001 )
// NOTE that (delta/epsilon_est) needs to be a positive integer!

// model problem information
struct ModelProblemData
  : public IModelProblemData
{
  static const bool has_exact_solution = true;

  ModelProblemData()
    : IModelProblemData(constants()) {
    if (constants().get("linear", true))
      DUNE_THROW(Dune::InvalidStateException, "problem eight is entirely nonlinear, but problem.linear was true");
  }

  //! \copydoc IModelProblemData::getMacroGridFile()
  inline std::string getMacroGridFile() const {
    return ("../dune/multiscale/grids/macro_grids/elliptic/unit_cube.dgf");
  }

  //! \copydoc IModelProblemData::getRefinementLevelReferenceProblem()
  inline int getRefinementLevelReferenceProblem() const {
    return 18;
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
    y = 0.0;
    y += 2.0 * ( x[0] + x[1] - pow(x[0], 2.0) - pow(x[1], 2.0) );
    y -= 12.0 * pow( (2.0 * x[0]) - 1.0, 2.0 ) * pow( (x[1] * x[1]) - x[1], 3.0 );
    y -= 12.0 * pow( (2.0 * x[1]) - 1.0, 2.0 ) * pow( (x[0] * x[0]) - x[0], 3.0 );
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

  double additivePart(const DomainType& x, const int i, const int j) const {
    double y = 0.0;

    y -= (x[0] + x[1]) * cos(2.0 * M_PI * x[i] / constants().epsilon) * sin(2.0 * M_PI * x[j] / constants().epsilon);

    double helper1 = 1.0;
    helper1 *= ( 2.0 + sin(2.0 * M_PI * (x[0] + x[1]) / constants().epsilon) );

    double helper2 = 1.0;
    helper2 *= 3.0;
    helper2 *= ( (2.0 * x[i] - 1.0) * (x[j] * x[j] - x[j]) )
               + ( (x[0] + x[1]) * cos(2.0 * M_PI * x[i] / constants().epsilon) * sin(2.0 * M_PI * x[j] / constants().epsilon) );
    helper2 *= (2.0 * x[i] - 1.0) * (x[j] * x[j] - x[j]) * (x[0] + x[1]) * cos(2.0 * M_PI * x[i] / constants().epsilon) * sin(
      2.0 * M_PI * x[j] / constants().epsilon);
    helper2 += pow( (x[0] + x[1]) * cos(2.0 * M_PI * x[i] / constants().epsilon) * sin(2.0 * M_PI * x[j] / constants().epsilon), 3.0 );

    helper1 *= helper2;
    y -= helper1;
    y -= sin(2.0 * M_PI * (x[0] + x[1]) / constants().epsilon) * pow( (2.0 * x[i]) - 1.0, 3.0 ) * pow( (x[j] * x[j]) - x[j], 3.0 );
    return y;
  } // additivePart

  // instantiate all possible cases of the evaluate-method:

  // (diffusive) flux = A^{\epsilon}( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& x,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const {
    double coefficient = 2.0 + sin(2.0 * M_PI * (x[0] + x[1]) / constants().epsilon);

    flux[0][0] = gradient[0][0] + ( coefficient * pow(gradient[0][0], 3.0) );
    flux[0][0] -= additivePart(x, 0, 1);
    flux[0][0] *= (-1.0);

    flux[0][1] = gradient[0][1] + ( coefficient * pow(gradient[0][1], 3.0) );
    flux[0][1] -= additivePart(x, 1, 0);
    flux[0][1] *= (-1.0);
  } // diffusiveFlux

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& x,
                             const JacobianRangeType& position_gradient,
                             const JacobianRangeType& direction_gradient,
                             JacobianRangeType& flux) const {
    double coefficient = 2.0 + sin(2.0 * M_PI * (x[0] + x[1]) / constants().epsilon);

    flux[0][0] = direction_gradient[0][0]
                 * ( 1.0 + 3.0 * coefficient * pow(position_gradient[0][0], 2.0) );
    flux[0][1] = direction_gradient[0][1]
                 * ( 1.0 + 3.0 * coefficient * pow(position_gradient[0][1], 2.0) );

    flux[0][0] *= (-1.0);
    flux[0][1] *= (-1.0);
  } // jacobianDiffusiveFlux

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate' method of the Diffusion class! See 'problem_specification.hh' for details.");
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
                     const JacobianRangeType& /*gradient*/,
                     JacobianRangeType& /*flux*/) const {
    DUNE_THROW(Dune::NotImplemented, "No homogenization available");
  }

  // the jacobian matrix (JA^{\epsilon}) of the diffusion operator A^{\epsilon} at the position "\nabla v" in direction
  // "nabla w", i.e.
  // jacobian diffusiv flux = JA^{\epsilon}(\nabla v) nabla w:

  // jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
  void jacobianDiffusiveFlux(const DomainType& /*x*/,
                             const JacobianRangeType& /*position_gradient*/,
                             const JacobianRangeType& /*direction_gradient*/,
                             JacobianRangeType& /*flux*/) const {
    if (DSC_CONFIG_GET("problem.linear", true))
      DUNE_THROW(Dune::NotImplemented, "Not yet implemented.");
    else
      DUNE_THROW(Dune::NotImplemented, "Nonlinear example not yet implemented.");
  } // jacobianDiffusiveFlux

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate' method of the Diffusion class! See 'problem_specification.hh' for details.");
  }
};

CONSTANTFUNCTION(MassTerm,  0.00001)
NULLFUNCTION(DefaultDummyFunction)

//! Exact solution (typically it is unknown)
template< class FunctionSpaceImp >
class ExactSolution
  : public Dune::Fem::Function< FunctionSpaceImp, ExactSolution< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef ExactSolution< FunctionSpaceType >                 ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;
  // essentially: 'DomainFieldType' is the type of an entry of a domain-element.
  // But: it is also used if 'u' (the exact solution) has a time-dependency ('u = u(x,t)').
  // This makes sense since the time-dependency is a one-dimensional element of the 'DomainType' and is therefor also an
  // entry of a domain-element.

public:
  ExactSolution(){}

  // in case 'u' has NO time-dependency use the following method:
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = (-1.0) * ( (x[0] * x[0]) - x[0] ) * ( (x[1] * x[1]) - x[1] );
  }

  // in case 'u' HAS a time-dependency use the following method:
  // unfortunately GRAPE requires both cases of the method 'evaluate' to be
  // instantiated
  inline void evaluate(const DomainType& x,
                       const TimeType& /*timedummy*/,
                       RangeType& y) const {
    evaluate(x, y);
  }

  inline void evaluateJacobian(const DomainType& , typename BaseType::JacobianRangeType& ) const {
    DUNE_THROW(Dune::NotImplemented, "Dummy body for all-problem compile");
  }
};
} // namespace Eight {
}

#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_EIGHT
