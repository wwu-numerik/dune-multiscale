#include "elliptic_homogenizer.hh"

namespace Dune {
namespace Multiscale {

Homogenizer::Homogenizer(const std::string &filename) : filename_(filename) {}

double Homogenizer::getEntry(const Homogenizer::TransformTensorType &tensor, const Homogenizer::PeriodicDiscreteFunctionSpaceType &periodicDiscreteFunctionSpace, const Homogenizer::PeriodicDiscreteFunctionType &w_i, const Homogenizer::PeriodicDiscreteFunctionType &w_j, const int &i, const int &j) const {
  double a_ij_hom = 0;

  for (const auto& entity : periodicDiscreteFunctionSpace) {
    auto localW_i = w_i.localFunction(entity);
    auto localW_j = w_j.localFunction(entity);

    // create quadrature for given geometry type
    //!\TODO soll das WIRLKLICH pold order 2 sein?
    const auto quadrature = make_quadrature(entity, periodicDiscreteFunctionSpace, 2);

    // get geoemetry of entity
    const auto& geometry = entity.geometry();

    // integrate
    for (const auto localQuadPoint : DSC::valueRange(quadrature.nop())) {
      RangeType localIntegral = 0;

      PeriodicJacobianRangeType grad_w_i;
      localW_i.jacobian(quadrature[localQuadPoint], grad_w_i);

      PeriodicJacobianRangeType grad_w_j;
      localW_j.jacobian(quadrature[localQuadPoint], grad_w_j);

      // local (barycentric) coordinates (with respect to cell grid entity)
      const auto& local_point = quadrature.point(localQuadPoint);

      // global point in the unit cell Y
      const auto global_point_in_Y = geometry.global(local_point);

      PeriodicJacobianRangeType direction_i;
      for (int k = 0; k < dimension; ++k) {
        direction_i[0][k] = grad_w_i[0][k];
        if (k == i) {
          direction_i[0][k] += 1.0;
        }
      }

      PeriodicJacobianRangeType direction_j;
      for (int k = 0; k < dimension; ++k) {
        direction_j[0][k] = 0.0;
        if (k == j) {
          direction_j[0][k] += 1.0;
        }
      }

      PeriodicJacobianRangeType flux_i;
      tensor.diffusiveFlux(global_point_in_Y, direction_i, flux_i);

      localIntegral = (flux_i[0] * direction_j[0]);

      const double entityVolume =
          quadrature.weight(localQuadPoint) * geometry.integrationElement(quadrature.point(localQuadPoint));

      a_ij_hom += entityVolume * localIntegral;
    }
  }

  return a_ij_hom;
}

Homogenizer::HomTensorType Homogenizer::getHomTensor(const Homogenizer::TensorType &tensor) const {
  HomTensorType a_hom;

  // to solve cell problems, we always need to use a perforated unit cube as domain:
  GridPtr<GridType> periodicgridptr(filename_);

  Dune::Fem::GlobalRefine::apply(*periodicgridptr, 10);

  PeriodicGridPartType periodicGridPart(*periodicgridptr);

  PeriodicDiscreteFunctionSpaceType periodicDiscreteFunctionSpace(periodicGridPart);

  // to avoid confusions:
  const DummySpaceType dummySpace(periodicGridPart);
  // (sometimes periodicDiscreteFunctionSpace is only a dummy)

  //! define the type of the corresponding solutions ( discrete functions of the type 'DiscreteFunctionType'):

  PeriodicDiscreteFunctionType cellSolution_0("cellSolution 0", periodicDiscreteFunctionSpace);
  cellSolution_0.clear();

  PeriodicDiscreteFunctionType cellSolution_1("cellSolution 1", periodicDiscreteFunctionSpace);
  cellSolution_1.clear();

  PeriodicDiscreteFunctionType rhs_0("rhs_0", periodicDiscreteFunctionSpace);
  rhs_0.clear();

  PeriodicDiscreteFunctionType rhs_1("rhs_1", periodicDiscreteFunctionSpace);
  rhs_1.clear();

  const RangeType lambda = 1e-07;
  // we need solve several cell problems with a periodic boundary condition. To fix the solution we need some kind of
  // additional condition. Therefor we solve
  // \lambda w - \div A \nabla w = rhs        instead of      - \div A \nabla w = rhs

  const TransformTensorType tensor_transformed(tensor);

  // if we have some additional source term (-div G), define:
  const CellSourceType G_0(periodicDiscreteFunctionSpace, tensor_transformed, 0); // 0'th cell problem
  const CellSourceType G_1(periodicDiscreteFunctionSpace, tensor_transformed, 1); // 1'th cell problem
  // - div ( A \nabla u^{\epsilon} ) = f - div G

  // quite a dummy. It's always f = 0
  const ZeroFunctionType zero;

  //! build the left hand side (lhs) of the problem
  const std::unique_ptr<const Problem::LowerOrderBase> mass(new MassWeightType(lambda));
  // define mass (just for cell problems \lambda w - \div A \nabla w = rhs)
  const EllipticOperatorType discrete_cell_elliptic_op(periodicDiscreteFunctionSpace, tensor_transformed, mass);

  LinearOperatorType lhsMatrix("Cell Problem Stiffness Matrix", periodicDiscreteFunctionSpace,
                               periodicDiscreteFunctionSpace);
  discrete_cell_elliptic_op.assemble_matrix(lhsMatrix, false /*no boundary treatment*/);

  //! build the right hand side (rhs) of the problem

  // the same right hand side for HM and FEM methods:
  typedef RightHandSideAssembler<PeriodicDiscreteFunctionType> RhsAssembler;

  // Alternativly it is possible to call the RightHandSideAssembler with a second source Term '- div G':
  // RightHandSideAssembler< DiscreteFunctionType > rhsassembler( tensor , G );
  RhsAssembler::assemble<2 * PeriodicDiscreteFunctionSpaceType::polynomialOrder>(zero, G_0, rhs_0);
  RhsAssembler::assemble<2 * PeriodicDiscreteFunctionSpaceType::polynomialOrder>(zero, G_1, rhs_1);
  // solve the linear systems (with Bi-CG):

  const InverseLinearOperatorType fembiCG(lhsMatrix, 1e-8, 1e-8, 20000,
                                          DSC_CONFIG_GET("global.cgsolver_verbose", false));

  fembiCG(rhs_0, cellSolution_0);
  fembiCG(rhs_1, cellSolution_1);

  a_hom[0][0] = getEntry(tensor_transformed, periodicDiscreteFunctionSpace, cellSolution_0, cellSolution_0, 0, 0);
  a_hom[1][1] = getEntry(tensor_transformed, periodicDiscreteFunctionSpace, cellSolution_1, cellSolution_1, 1, 1);
  a_hom[0][1] = getEntry(tensor_transformed, periodicDiscreteFunctionSpace, cellSolution_0, cellSolution_1, 0, 1);
  a_hom[1][0] = a_hom[0][1];

  DSC_LOG_DEBUG << "A_homogenized[0][0] = " << a_hom[0][0] << std::endl << "A_homogenized[0][1] = " << a_hom[0][1]
      << std::endl << "A_homogenized[1][0] = " << a_hom[1][0] << std::endl
                                                              << "A_homogenized[1][1] = " << a_hom[1][1] << std::endl;

  return a_hom;
}

} // namespace Multiscale {
} // end namespace

