// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNEMS_HMM_CELL_SOLVER_HH
#define DUNEMS_HMM_CELL_SOLVER_HH

#include <dune/multiscale/hmm/hmm_traits.hh>
#include <dune/multiscale/common/traits.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/profiler.hh>
#include <dune/stuff/common/math.hh>

#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/solver/oemsolver.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

#include <dune/multiscale/tools/discretefunctionwriter.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>

//! --------------------- the essential cell problem solver class ----------------------------------
namespace Dune {
namespace Multiscale {
namespace HMM {

class CellProblemSolver
{
private:
  typedef typename HMMTraits::PeriodicDiscreteFunctionType PeriodicDiscreteFunctionImp;
  typedef typename CommonTraits::DiffusionType DiffusionImp;
  //! type of discrete functions
  typedef PeriodicDiscreteFunctionImp PeriodicDiscreteFunctionType;

  //! type of discrete function space
  typedef typename PeriodicDiscreteFunctionType::DiscreteFunctionSpaceType
  PeriodicDiscreteFunctionSpaceType;

  //! type of grid partition
  typedef typename PeriodicDiscreteFunctionSpaceType::GridPartType PeriodicGridPartType;

  //! type of grid
  typedef typename PeriodicDiscreteFunctionSpaceType::GridType PeriodicGridType;

  //! type of range vectors
  typedef typename PeriodicDiscreteFunctionSpaceType::RangeType RangeType;

  //! type of range vectors
  typedef typename PeriodicDiscreteFunctionSpaceType::DomainType DomainType;

  //! polynomial order of base functions
  enum { polynomialOrder = PeriodicDiscreteFunctionSpaceType::polynomialOrder };

  //! type of the (possibly non-linear) diffusion operator
  typedef DiffusionImp DiffusionType;

  struct CellMatrixTraits
  {
    typedef PeriodicDiscreteFunctionSpaceType                          RowSpaceType;
    typedef PeriodicDiscreteFunctionSpaceType                          ColumnSpaceType;
    typedef LagrangeMatrixSetup< false >     StencilType;
    typedef ParallelScalarProduct< PeriodicDiscreteFunctionSpaceType > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      //!TODO this works only with ISTL present
      typedef Dune::LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };
public:
  typedef SparseRowMatrixOperator< PeriodicDiscreteFunctionType, PeriodicDiscreteFunctionType,
                                   CellMatrixTraits > CellFEMMatrix;
private:
  //! OEMGMRESOp //OEMBICGSQOp // OEMBICGSTABOp
  typedef OEMBICGSTABOp< PeriodicDiscreteFunctionType, CellFEMMatrix > InverseCellFEMMatrix;

private:
  const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace_; // Referenz &, wenn & verwendet, dann unten:
  const DiffusionType& diffusion_;
  static const std::string subdir_; /*"cell_problems"*/

public:

  //! constructor - with diffusion operator A^{\epsilon}(x)
  CellProblemSolver(const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                    const DiffusionType& diffusion_operator)
    : periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
      , diffusion_(diffusion_operator)
  {}

private:
  //! ----------- method: solve cell problem ------------------------------------------
  void solvecellproblem(const typename CommonTraits::DiscreteFunctionType::DiscreteFunctionSpaceType::JacobianRangeType& gradient_PHI_H,
                        // the barycenter x_T of a macro grid element 'T'
                        const DomainType& globalQuadPoint,
                        PeriodicDiscreteFunctionType& cell_problem_solution) const;

  /** -------- method: solve the jacobian corrector cell problem ----------------------------
   * Problem to determine the Jacobian operator of the correction operator
   * ( the correction operator Q is defined via cell problems, here we determine the corresponding derivative DQ )
   * (see paper what it means and where it comes from. Note that it is only required for the nonlinear case)
   * \arg gradient_PHI_H gradient of macroscopic base function
   * \arg gradient of the macroscopic function from the last iteration step
   * \arg gradient_y of the corrector of the macroscopic function from the last iteration step
   * \arg the barycenter x_T of a macro grid element 'T'
   **/
  void solve_jacobiancorrector_cellproblem(
    const typename CommonTraits::DiscreteFunctionType::DiscreteFunctionSpaceType::JacobianRangeType& gradient_PHI_H,
    const typename CommonTraits::DiscreteFunctionType::DiscreteFunctionSpaceType::JacobianRangeType& grad_old_coarse_function,
    const PeriodicDiscreteFunctionType& corrector_of_old_coarse_function,
    const DomainType& globalQuadPoint,
    PeriodicDiscreteFunctionType& jac_cor_cell_problem_solution) const;

public:
  /** two methods for solving and saving the solutions of the cell problems.
   * 1. save solutions for the whole set of macroscopic base function
   * (in general for the case of a linear diffusion operator)
   * 2. save the solutions for a fixed discrete function
   * (in general for the case of a nonlinear diffusion operator)
   *! ---- method: solve and save the cell problems for the set of macroscopic base functions -----
   * here we need a 'cell problem numbering manager' to determine the number of the cell problem
   * (a combination of number of entity and number of local base function)
   * Structure:
   * Struktur der Indizierung fuer das Abspeichern der Loesungen der Zellprobleme:
   * wir loesen Zellprobleme fuer jede Entity des Makro-Grids und jede Basisfunktion, die einen nichtleeren support auf
   * dieser Entity besitzt, also schematisch:
   * Sei n=0,...,N einer Durchnumerierung der Entitys und i=0,...I_n eine zu einer festen Entity gehoerende Nummerierung
   * der Basisfunktionen mit nicht-leeren support.
   * Die Durchnummerierung der Loesungen der Zellprobleme k=0,...,K ist dann gegeben durch: k(n,i_n) = ( sum_(l=0)^(n-1)
   * ( I_l + 1) ) + i_n
   * NOTE: es verhaelt sich NICHT wie die vorhandene Methode mapToGlobal(entity,i) ! (die gibt die globale Nummer der
   * Basisfunktion zurueck, es gibt aber  deutlich mehr Zellprobleme zum Loesen!
   * (das wird aber alles im Hintergrund vom 'cell problem numbering manager')

   * compute and save solutions of the cell problems for the base function set of the 'discreteFunctionSpace'
   * requires cell problem numbering manager
   **/
  void saveTheSolutions_baseSet(
    const typename CommonTraits::DiscreteFunctionType::DiscreteFunctionSpaceType& discreteFunctionSpace,
    const CellProblemNumberingManager& cp_num_manager   // just to check, if we use the correct numeration
      ) const;


  /** this method does not require a cell problem numbering manager
   * (it uses the standard counter for entities, provided by the corrsponding iterator)
   * compute and save solutions of the cell problems for a fixed macroscopic discrete function
   * (in gerneral it is the macro solution from the last iteration step)
   **/
  void saveTheSolutions_discFunc(
    const CommonTraits::DiscreteFunctionType& macro_discrete_function) const;

  //! compute and save solutions of the jacobian corrector cell problems for the base function set of the
  //! 'discreteFunctionSpace' and for a fixed macroscopic discrete function
  //! requires cell problem numbering manager
  void saveTheJacCorSolutions_baseSet_discFunc(
    const CommonTraits::DiscreteFunctionType& macro_discrete_function,
    const CellProblemNumberingManager& cp_num_manager   // just to check, if we use the correct numeration
    ) const;
}; // end class




} //namespace HMM {
} //namespace Multiscale {
} //namespace Dune {


#endif // DUNEMS_HMM_CELL_SOLVER_HH
