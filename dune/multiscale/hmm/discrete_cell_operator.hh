// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DiscreteEllipticCellProblem_HH
#define DiscreteEllipticCellProblem_HH

#ifdef HAVE_CMAKE_CONFIG
 #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
 #include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <dune/common/fmatrix.hh>
#include <boost/noncopyable.hpp>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>

#include <dune/fem/operator/2order/lagrangematrixsetup.hh>

#include <dune/multiscale/hmm/cell_problem_solver.hh>
#include <dune/multiscale/hmm/hmm_traits.hh>
#include <dune/multiscale/problems/elliptic/selector.hh>

namespace Dune {
namespace Multiscale {
namespace HMM {


class DiscreteCellProblemOperator
  : public Operator< typename HMMTraits::PeriodicDiscreteFunctionType::RangeFieldType,
                     typename HMMTraits::PeriodicDiscreteFunctionType::RangeFieldType, typename HMMTraits::PeriodicDiscreteFunctionType,
                     typename HMMTraits::PeriodicDiscreteFunctionType >
    , boost::noncopyable
{
  typedef typename HMMTraits::PeriodicDiscreteFunctionType PeriodicDiscreteFunctionImp;
  typedef typename Problem::Diffusion DiffusionImp;

  typedef PeriodicDiscreteFunctionImp DiscreteFunction;
  typedef DiffusionImp                DiffusionModel;

  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpace;

  typedef typename DiscreteFunctionSpace::GridPartType   GridPart;
  typedef typename DiscreteFunctionSpace::GridType       GridType;
  typedef typename DiscreteFunctionSpace::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpace::DomainType DomainType;
  typedef typename DiscreteFunctionSpace::RangeType  RangeType;
  typedef typename DiscreteFunctionSpace::JacobianRangeType
  JacobianRangeType;

  static const int dimension = GridPart::GridType::dimension;
  static const int polynomialOrder = DiscreteFunctionSpace::polynomialOrder;

  typedef typename DiscreteFunction::LocalFunctionType LocalFunction;

  typedef typename DiscreteFunctionSpace::BaseFunctionSetType                   BaseFunctionSet;
  typedef typename DiscreteFunctionSpace::LagrangePointSetType                  LagrangePointSet;
  typedef typename LagrangePointSet::Codim< 1 >::SubEntityIteratorType FaceDofIterator;

  typedef typename DiscreteFunctionSpace::IteratorType Iterator;
  typedef typename Iterator::Entity                    Entity;
  typedef typename Entity::Geometry                    Geometry;

  typedef typename GridPart::IntersectionIteratorType IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;

  typedef CachingQuadrature< GridPart, 0 > Quadrature;

public:
  DiscreteCellProblemOperator(const DiscreteFunctionSpace& periodicDiscreteFunctionSpace,
                              const DiffusionModel& diffusion_op)
    : periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
      , diffusion_operator_(diffusion_op)
  {}

public:
  /**
   * @brief dummy implementation of "operator()"
   * @param w effect of the discrete operator on 'u'
   */
  virtual void operator()(const DiscreteFunction& u, DiscreteFunction& w) const;

  /**
   ! stiffness matrix for a linear elliptic diffusion operator
   we obtain entries of the following kind
   (cell problem for the macro grid element 'T' and for the base-function '\Phi_H',
   x_T denotes the barycenter of T, \delta denotes the cell size )
   \int_Y A_h^{\eps}(t,x_T + \delta*y) \nabla phi_h_i(y) \cdot \nabla phi_h_j(y)
   + CELL_MASS_WEIGHT * \int_Y phi_h_i(y) \phi_h_j(y)
   (the second summand yields an artificical mass term to guarantee uniqueness and existence
   for the problem with periodic boundary condition.
   This is an alternative to the 'average zero' condition.)
  **/
  void assemble_matrix(const DomainType& x_T, typename CellProblemSolver::CellFEMMatrix& global_matrix) const;

  /**
   *  ! stiffness matrix for a non linear elliptic diffusion operator

 here we obtain the contribution of the jacobian matrix of the diffusion operator evaluated in the gradient of a
 certain discrete function (in case of the Newton method, it is the preceeding iterate v_h^{(n-1)} )

 Let PHI_H denote a macrocopic discrete function (that we want to reconstruct)
 we obtain entries of the following kind:
 (cell problem for the macro grid element 'T')

 \int_Y JA^{\eps}( x_T + \delta*y, \nabla_x PHI_H(x_T) + \nabla_y v_h^{(n-1)} ) \nabla phi_h_i(y) \cdot \nabla
 phi_h_j(y)
 + CELL_MASS_WEIGHT * \int_Y phi_h_i(y) \phi_h_j(y)

 (here, JA^{\eps} denotes the jacobian matrix of the diffusion operator A^{\eps},
 x_T denotes the barycenter of T, \delta denotes the cell size )
   * \param x_T macroscopic quadrature point
   * \param old_fine_function the microscopic function (fine-scale correction) from the last iteration step of the Newton method
   * \param grad_coarse_function the gradient of the macroscopic function (that we want to reconstruct) evaluated in x_T
   */
  void assemble_jacobian_matrix(const DomainType& x_T,
                                const JacobianRangeType& grad_coarse_function,
                                const DiscreteFunction& old_fine_function,
                                typename CellProblemSolver::CellFEMMatrix& global_matrix) const;

  // begin group the "right hand side assembler methods"
  /**
   * assemble method for the case of a linear diffusion operator
   (in this case, no Newton method is required, which is why there is no dependency on an old fine-scale discrete
   function / old iteration step )

   we compute the following entries for each fine-scale base function phi_h_i:
   - \int_Y A^{\eps}( x_T + \delta*y ) \nabla_x PHI_H(x_T) \cdot \nabla_y phi_h_i(y)
   * @brief assembleCellRHS_linear
   * @param x_T the global quadrature point in the macro grid element T
   * @param grad_coarse_function \nabla_x \Phi_H(x_T) (the coarse function to reconstruct)
   * @param cell_problem_RHS rhs cell problem
   */
  void assembleCellRHS_linear( const DomainType& x_T,
    const JacobianRangeType& grad_coarse_function,
    DiscreteFunction& cell_problem_RHS) const;

  /**
   * assemble method for the case of a nonlinear diffusion operator
   (in this case a Newton method is required, which is why there is an additional dependency on an old fine-scale
   discrete function v_h^{(n-1)} / old iteration step )

   we compute the following entries for each fine-scale base function phi_h_i:
   - \int_Y A^{\eps}( x_T + \delta*y , \nabla_x PHI_H(x_T) + \nabla_y \nabla_y v_h^{(n-1)} ) \cdot \nabla_y phi_h_i(y)
   * @brief assembleCellRHS_nonlinear
   * @param x_T the global quadrature point in the macro grid element T
   * @param grad_coarse_function \nabla_x \Phi_H(x_T) :
      in the linear setting we typically reconstruct macroscopic base functions, in the non-linear setting we are in a
      more general setting.
      gradient of the coarse function, that we want to reconstruct:
   * @param old_fine_function old solution from the last iteration step
   * @param cell_problem_RHS rhs cell problem
   */
  void assembleCellRHS_nonlinear( const DomainType& x_T,
    const JacobianRangeType& grad_coarse_function,
    const DiscreteFunction& old_fine_function,
    DiscreteFunction& cell_problem_RHS) const;

  /**
   * @brief assemble_jacobian_corrector_cell_prob_RHS assemble method for the right hand side of the jacobian corrector cell problem
   * @param x_T the global quadrature point in the macro grid element T
   * @param grad_old_coarse_function gradient of the old coarse function (old means last iteration step)
   * @param corrector_of_old_coarse_function gradient of the corrector of the old coarse function
   * @param grad_coarse_base_function gradient of the current macroscopic base function
   * @param jac_corrector_cell_problem_RHS rhs cell problem
   */
  void assemble_jacobian_corrector_cell_prob_RHS( const DomainType& x_T,
    const JacobianRangeType& grad_old_coarse_function,
    const DiscreteFunction& corrector_of_old_coarse_function,
    const JacobianRangeType& grad_coarse_base_function,
    DiscreteFunction& jac_corrector_cell_problem_RHS) const;

  void printCellRHS(const DiscreteFunction& rhs) const;

  double normRHS(const DiscreteFunction& rhs) const;

private:
  const DiscreteFunctionSpace& periodicDiscreteFunctionSpace_;
  const DiffusionModel& diffusion_operator_;
};

} //namespace HMM {
} //namespace Multiscale {
} //namespace Dune {


#endif // #ifndef DiscreteElliptic_HH
