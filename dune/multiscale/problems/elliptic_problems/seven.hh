#ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SEVEN
#define DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SEVEN

#include <dune/fem/function/common/function.hh>
#include <dune/multiscale/problems/constants.hh>
#include <dune/multiscale/problems/base.hh>

//! ------------ Elliptic Problem 7 -------------------

// linear elliptic model problem - periodic setting

//! For more further details about the implementation see '../base.hh'
//! For details on the classes, see 'example.hh'

// Note that in the following, 'Imp' abbreviates 'Implementation'
namespace Problem {
namespace Seven {
// default value for epsilon (if not sprecified in the parameter file)
CONSTANTSFUNCTION( 0.01 )

// model problem information
struct ModelProblemData
  : public IModelProblemData
{
  static const bool has_exact_solution = false;

  ModelProblemData()
    : IModelProblemData(constants()) {
      if (!constants().get("linear", true))
        DUNE_THROW(Dune::InvalidStateException, "problem seven is entirely linear, but problem.linear was false");
      assert( constants_.epsilon != 0.0);
      if (constants().get("stochastic_pertubation", false) && !(this->problemAllowsStochastics()) )
         DUNE_THROW(Dune::InvalidStateException, "The problem does not allow stochastic perturbations. Please, switch the key off.");      
    }

  //! \copydoc IModelProblemData::getMacroGridFile()
  inline std::string getMacroGridFile() const {
    return("../dune/multiscale/grids/macro_grids/elliptic/cube_three.dgf");
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

//! ----------------- Definition of ' f ' ------------------------
CONSTANTFUNCTION(FirstSource, 1.0)
//! ----------------- End Definition of ' f ' ------------------------


//! ----------------- Definition of ' G ' ------------------------
NULLFUNCTION(SecondSource)
//! ----------------- End Definition of ' G ' ------------------------


//! ----------------- Definition of ' A ' ------------------------
// the linear diffusion operator A^{\epsilon}(x,\xi)=A^{\epsilon}(x) \xi
// A^{\epsilon} : \Omega × R² -> R²
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

   // (diffusive) flux = A^{\epsilon}( x , direction )
   // (typically direction is some 'gradient_of_a_function')

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

  template < class... Args >
  void evaluate( Args... ) const
  {
    DUNE_THROW(Dune::NotImplemented, "Inadmissible call for 'evaluate'");
  }

};
//! ----------------- End Definition of ' A ' ------------------------


//! ----------------- Definition of ' m ' ----------------------------
CONSTANTFUNCTION(MassTerm,  0.0)
//! ----------------- End Definition of ' m ' ------------------------


//! ----------------- Definition of some dummy -----------------------
NULLFUNCTION(DefaultDummyFunction)
//! ----------------- End Definition of some dummy -------------------


// Exact solution is unknown:
//! ----------------- Definition of ' u ' ----------------------------
NULLFUNCTION(ExactSolution)
//! ----------------- End Definition of ' u ' ------------------------


} // namespace Seven {
}


#endif // ifndef DUNE_ELLIPTIC_MODEL_PROBLEM_SPECIFICATION_HH_SEVEN
