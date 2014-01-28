// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_ERRORESTIMATER_HH
#define DUNE_ERRORESTIMATER_HH


#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/stuff/common/parameter/configcontainer.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/entity.hh>

//!NOTE: 'ErrorEstimator' requires an access to the 'ModelProblemData' class (typically defined in
// problem_specification.hh), which provides us with infomration about epsilon, delta, etc.

//! ------------------------- nonlinear elliptic ------------------------------------------

namespace Dune {
namespace Multiscale {
namespace HMM {

//! NOTE: das zweite Argument DiscreteFunctionImp sobald wie moeglich wieder loeschen. Alle 'Ableitungen' und Variablen
//! von DiscreteFunctionImp ebenfalls loeschen und die zugehörigen von PeriodicDiscreteFunctionImp ersetzen. Aktuell ist
//! es nur da, weil verschiedene Dinge wie der IntersectionIterator für das periodic Gridpart noch nicht implementiert
//! sind.
template <class PeriodicDiscreteFunctionImp, class DiscreteFunctionImp, class DiffusionImp>
class ErrorEstimator {
private:
  typedef DiffusionImp DiffusionOperatorType;

  //! Necessary typedefs for the PeriodicDiscreteFunctionImp:

  typedef PeriodicDiscreteFunctionImp PeriodicDiscreteFunctionType;

  typedef typename PeriodicDiscreteFunctionType::LocalFunctionType PeriodicLocalFunctionType;
  typedef typename PeriodicDiscreteFunctionType::DiscreteFunctionSpaceType PeriodicDiscreteFunctionSpaceType;
  typedef typename PeriodicDiscreteFunctionType::RangeType RangeType;
  typedef typename PeriodicDiscreteFunctionType::RangeFieldType TimeType;
  typedef typename PeriodicDiscreteFunctionType::DomainType DomainType;
  typedef typename PeriodicDiscreteFunctionSpaceType::JacobianRangeType PeriodicJacobianRangeType;

  typedef typename PeriodicDiscreteFunctionType::GridType::template Codim<0>::Entity PeriodicEntityType;
  typedef typename PeriodicDiscreteFunctionType::GridType::template Codim<0>::EntityPointer PeriodicEntityPointerType;
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::GridType::template Codim<0>::Entity EntityType;
  typedef typename DiscreteFunctionSpaceType::GridType::template Codim<0>::EntityPointer EntityPointerType;

  static const int dimension = PeriodicDiscreteFunctionType::GridType::dimension;
  static const int spacePolOrd = PeriodicDiscreteFunctionSpaceType::CommonTraits::polynomial_order;
  static const int maxnumOfBaseFct = 100;

private:
  const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace_;
  const DiscreteFunctionSpaceType& discreteFunctionSpace_;
  const DiscreteFunctionSpaceType& auxiliaryDiscreteFunctionSpace_;
  // an auxiliaryDiscreteFunctionSpace to get an Intersection Iterator for the periodicDiscreteFunctionSpace
  // (for the periodic grid partition there is no usable intersection iterator implemented, therefore we use the
  // intersection iterator for the corresponding non-periodic grid partition (this not efficient and increases the
  // estimated error, but it works)

  const DiffusionOperatorType& diffusion_;

public:
  ErrorEstimator(const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                 const DiscreteFunctionSpaceType& discreteFunctionSpace,
                 const DiscreteFunctionSpaceType& auxiliaryDiscreteFunctionSpace,
                 const DiffusionOperatorType& diffusion)
    : periodicDiscreteFunctionSpace_(periodicDiscreteFunctionSpace)
    , discreteFunctionSpace_(discreteFunctionSpace)
    , auxiliaryDiscreteFunctionSpace_(auxiliaryDiscreteFunctionSpace)
    , diffusion_(diffusion) {}

  RangeType getH(const EntityType& entity) const { return DSG::entityDiameter(entity); } // getH

  // return:  H_T^4 ||f||_{L^2(T)}^2
  template <class SourceType>
  RangeType indicator_f(const SourceType& f, const EntityType& entity) const {
    const auto entityQuadrature = DSFe::make_quadrature(entity, periodicDiscreteFunctionSpace_);

    // get geoemetry of entity
    const auto& geometry = entity.geometry();

    RangeType H_T = getH(entity);

    RangeType y(0);
    RangeType local_indicator(0);

    const auto quadratureNop = entityQuadrature.nop();
    for (size_t quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint) {
      const double weight = entityQuadrature.weight(quadraturePoint) *
                            geometry.integrationElement(entityQuadrature.point(quadraturePoint));

      f.evaluate(geometry.global(entityQuadrature.point(quadraturePoint)), y);
      y = y * y;

      local_indicator += weight * y;
    }

    local_indicator *= pow(H_T, 4.0);

    return local_indicator;
  } // indicator_f

  // \eta_T^{app}
  RangeType indicator_app_1(const EntityType& entity, const DiscreteFunctionType& u_H,
                            const PeriodicDiscreteFunctionType& corrector_u_H_on_entity) const {
    RangeType local_indicator(0.0);

    const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);
    const double epsilon_estimated = DSC_CONFIG_GET("hmm.epsilon_guess", 1.0f);

    const auto& globalEntityGeometry = entity.geometry();
    const auto& x_T = globalEntityGeometry.center();
    const auto& x_T_local = globalEntityGeometry.local(x_T);
    const auto entityVolume = globalEntityGeometry.volume();

    // int_Y A_h^{\epsilon} - A^{\epsilob}
    RangeType difference[dimension];

    for (int k = 0; k < dimension; ++k)
      difference[k] = 0.0;

    // \nabla u_H(x_T)
    const auto u_H_local = u_H.localFunction(entity);
    JacobianRangeType gradient_u_H(0.);
    u_H_local.jacobian(x_T_local, gradient_u_H);

    for (const auto& micro_entity : periodicDiscreteFunctionSpace_) {
      // one quadrature formula ( A_h^{\eps}(y)=A^{\eps}(y_s) vs. A^{\eps}(y) )
      const auto center_local = micro_entity.geometry().local(micro_entity.geometry().center());
      const auto high_order_quadrature = DSFe::make_quadrature(micro_entity, periodicDiscreteFunctionSpace_);

      // Q_h(u_H)(x_T,y) on the micro entity:
      auto loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(micro_entity);
      JacobianRangeType gradient_Q_u_H_x_T(0.);
      loc_Q_u_H_x_T.jacobian(center_local, gradient_Q_u_H_x_T);

      // S denotes the micro grid element (i.e. 'micro_entity')
      const auto& geometry_S = micro_entity.geometry();
      // y_S denotes the barycenter of the micro grid element S:
      const auto& y_S = geometry_S.center();

      // to evaluate A^{\epsilon}_h:
      DomainType globalPoint_center;
      JacobianRangeType direction_of_diffusion;
      for (int k = 0; k < dimension; ++k) {
        globalPoint_center[k] = x_T[k] + (delta * y_S[k]);
        direction_of_diffusion[0][k] = gradient_u_H[0][k] + gradient_Q_u_H_x_T[0][k];
      }

      JacobianRangeType diffusive_flux_y_S;
      diffusion_.diffusiveFlux(globalPoint_center, direction_of_diffusion, diffusive_flux_y_S);

      double cutting_function = 1.0;
      for (int k = 0; k < dimension; ++k) {
        // is the current quadrature point in the relevant cell?
        // Y = [ -0.5 , 0.5 ]^dimension
        if (fabs(y_S[k]) > (0.5 * (epsilon_estimated / delta))) {
          cutting_function *= 0.0;
        }
      }

      const auto numQuadraturePoints = high_order_quadrature.nop();
      for (size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint) {
        // local (barycentric) coordinates (with respect to entity)
        const auto& y_barycentric = high_order_quadrature.point(microQuadraturePoint);

        // current quadrature point y in global coordinates
        auto y = geometry_S.global(y_barycentric);

        const double weight_micro_quadrature =
            high_order_quadrature.weight(microQuadraturePoint) * geometry_S.integrationElement(y_barycentric);

        // to evaluate A^{\epsilon}:
        DomainType globalPoint;
        for (int k = 0; k < dimension; ++k)
          globalPoint[k] = (delta * y[k]) + x_T[k];

        // diffusive flux in y (in direction \nabla u_H(x_T) + \nabla_y Q_h^T(u_H)(y_S) )
        JacobianRangeType diffusive_flux_y;
        diffusion_.diffusiveFlux(globalPoint, direction_of_diffusion, diffusive_flux_y);

        for (int k = 0; k < dimension; ++k)
          difference[k] +=
              cutting_function * weight_micro_quadrature * (diffusive_flux_y[0][k] - diffusive_flux_y_S[0][k]);
      }
    }

    for (int k = 0; k < dimension; ++k)
      local_indicator += pow(difference[k], 2.0) * entityVolume;

    return local_indicator;
  } // end of method

  // \bar{\eta}_T^{app}
  RangeType indicator_app_2(const EntityType& entity, const DiscreteFunctionType& u_H,
                            const PeriodicDiscreteFunctionType& corrector_u_H_on_entity) const {
    RangeType local_indicator(0.0);

    const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);
    const auto& globalEntityGeometry = entity.geometry();
    const auto& x_T = globalEntityGeometry.center();
    const auto& x_T_local = globalEntityGeometry.local(x_T);
    const auto entityVolume = globalEntityGeometry.volume();

    // int_Y A_h^{\epsilon} - A^{\epsilob}
    RangeType difference[dimension];

    for (int k = 0; k < dimension; ++k)
      difference[k] = 0.0;

    // \nabla u_H(x_T)
    const auto u_H_local = u_H.localFunction(entity);
    JacobianRangeType gradient_u_H(0.);
    u_H_local.jacobian(x_T_local, gradient_u_H);

    // iterator over the elements of the periodic micro grid:
    for (const auto& micro_entity : periodicDiscreteFunctionSpace_) {
      // two quadrature formulas ( A_h^{\eps}(y)=A^{\eps}(y_s) vs. A^{\eps}(y) )
      const auto high_order_quadrature = DSFe::make_quadrature(micro_entity, periodicDiscreteFunctionSpace_);

      // Q_h(u_H)(x_T,y) on the micro entity:
      const auto loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(micro_entity);
      JacobianRangeType gradient_Q_u_H_x_T(0.);
      loc_Q_u_H_x_T.jacobian(x_T_local, gradient_Q_u_H_x_T);

      // S denotes the micro grid element (i.e. 'micro_entity')
      const auto& geometry_S = micro_entity.geometry();
      // y_S denotes the barycenter of the micro grid element S:
      const DomainType& y_S = geometry_S.center();

      // to evaluate A^{\epsilon}_h:
      DomainType globalPoint_center;
      JacobianRangeType direction_of_diffusion;
      for (int k = 0; k < dimension; ++k) {
        globalPoint_center[k] = x_T[k] + (delta * y_S[k]);
        direction_of_diffusion[0][k] = gradient_u_H[0][k] + gradient_Q_u_H_x_T[0][k];
      }

      JacobianRangeType diffusive_flux_y_S;
      diffusion_.diffusiveFlux(globalPoint_center, direction_of_diffusion, diffusive_flux_y_S);

      const auto numQuadraturePoints = high_order_quadrature.nop();
      for (size_t microQuadraturePoint = 0; microQuadraturePoint < numQuadraturePoints; ++microQuadraturePoint) {
        // local (barycentric) coordinates (with respect to entity)
        const auto& y_barycentric = high_order_quadrature.point(microQuadraturePoint);
        // current quadrature point y in global coordinates
        const auto y = geometry_S.global(y_barycentric);
        const double weight_micro_quadrature =
            high_order_quadrature.weight(microQuadraturePoint) * geometry_S.integrationElement(y_barycentric);

        // to evaluate A^{\epsilon}:
        DomainType globalPoint;
        for (int k = 0; k < dimension; ++k)
          globalPoint[k] = (delta * y[k]) + x_T[k];

        // diffusive flux in y (in direction \nabla u_H(x_T) + \nabla_y Q_h^T(u_H)(y_S) )
        JacobianRangeType diffusive_flux_y;
        diffusion_.diffusiveFlux(globalPoint, direction_of_diffusion, diffusive_flux_y);

        for (int k = 0; k < dimension; ++k)
          difference[k] += weight_micro_quadrature * (diffusive_flux_y[0][k] - diffusive_flux_y_S[0][k]);
      }
    }

    for (int k = 0; k < dimension; ++k)
      local_indicator += pow(difference[k], 2.0) * entityVolume;

    return local_indicator;
  } // end of method

  // (1/2)*H_E^3*\eta_E^{res}
  // (1/2) since we visit every entity twice (inner and outer)
  template <class IntersectionType>
  RangeType indicator_res_E(const IntersectionType& intersection, const DiscreteFunctionType& u_H,
                            const PeriodicDiscreteFunctionType& corrector_u_H_on_inner_entity,
                            const PeriodicDiscreteFunctionType& corrector_u_H_on_outer_entity) const {
    RangeType local_indicator(0.0);
    const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);
    const double epsilon_estimated = DSC_CONFIG_GET("hmm.epsilon_guess", 1.0f);

    auto inner_it = intersection.inside();
    auto outer_it = intersection.outside();

    const auto& inner_entity = *inner_it;
    const auto& outer_entity = *outer_it;

    const auto unitOuterNormal = intersection.centerUnitOuterNormal();

    const auto& faceGeometry = intersection.geometry();
    // H_E (= |E|):
    const auto edge_length = faceGeometry.volume();

    // jump = innerValue - outerValue
    RangeType innerValue = 0.0;
    RangeType outerValue = 0.0;

    const auto& globalInnerEntityGeometry = inner_entity.geometry();
    const auto& globalOuterEntityGeometry = outer_entity.geometry();
    const auto& x_T_inner = globalInnerEntityGeometry.center();
    const auto& x_T_outer = globalOuterEntityGeometry.center();
    const auto& x_T_inner_local = globalInnerEntityGeometry.local(x_T_inner);
    const auto& x_T_outer_local = globalOuterEntityGeometry.local(x_T_outer);

    // \nabla u_H(x_T) (on the inner element T)
    auto inner_u_H_local = u_H.localFunction(inner_entity);
    JacobianRangeType gradient_inner_u_H(0.);

    inner_u_H_local.jacobian(x_T_inner_local, gradient_inner_u_H);

    // \nabla u_H(x_T) (on the outer element \bar{T})
    auto outer_u_H_local = u_H.localFunction(outer_entity);
    JacobianRangeType gradient_outer_u_H(0.);
    outer_u_H_local.jacobian(x_T_outer_local, gradient_outer_u_H);

    for (const auto& micro_entity : periodicDiscreteFunctionSpace_) {
      const auto center_local = micro_entity.geometry().local(micro_entity.geometry().center());
      // Q_h(u_H)(x_T,y) on the micro entity:
      auto loc_Q_u_H_x_T_inner = corrector_u_H_on_inner_entity.localFunction(micro_entity);
      JacobianRangeType gradient_Q_u_H_x_T_inner(0.);
      loc_Q_u_H_x_T_inner.jacobian(center_local, gradient_Q_u_H_x_T_inner);

      // Q_h(u_H)(x_{bar{T}},y) on the micro entity:
      auto loc_Q_u_H_x_T_outer = corrector_u_H_on_outer_entity.localFunction(micro_entity);
      JacobianRangeType gradient_Q_u_H_x_T_outer(0.);
      loc_Q_u_H_x_T_outer.jacobian(center_local, gradient_Q_u_H_x_T_outer);

      // S denotes the micro grid element (i.e. 'micro_entity')
      const auto& geometry_S = micro_entity.geometry();

      // y_S denotes the barycenter of the micro grid element S:
      const auto& y_S = geometry_S.center();

      // to evaluate A^{\epsilon}_h:
      DomainType globalPoint_inner;
      JacobianRangeType direction_of_diffusion_inner;
      for (int k = 0; k < dimension; ++k) {
        globalPoint_inner[k] = x_T_inner[k] + (delta * y_S[k]);
        direction_of_diffusion_inner[0][k] = gradient_inner_u_H[0][k] + gradient_Q_u_H_x_T_inner[0][k];
      }

      // to evaluate A^{\epsilon}_h:
      DomainType globalPoint_outer;
      JacobianRangeType direction_of_diffusion_outer;
      for (int k = 0; k < dimension; ++k) {
        globalPoint_outer[k] = x_T_outer[k] + (delta * y_S[k]);
        direction_of_diffusion_outer[0][k] = gradient_outer_u_H[0][k] + gradient_Q_u_H_x_T_outer[0][k];
      }

      JacobianRangeType diffusive_flux_y_S_inner;
      diffusion_.diffusiveFlux(globalPoint_inner, direction_of_diffusion_inner, diffusive_flux_y_S_inner);

      JacobianRangeType diffusive_flux_y_S_outer;
      diffusion_.diffusiveFlux(globalPoint_outer, direction_of_diffusion_outer, diffusive_flux_y_S_outer);

      double cutting_function = 1.0;
      for (int k = 0; k < dimension; ++k) {
        // is the current quadrature point in the relevant cell?
        // Y = [ -0.5 , 0.5 ]^dimension
        if (fabs(y_S[k]) > (0.5 * (epsilon_estimated / delta))) {
          cutting_function *= 0.0;
        }
      }

      const double weight_micro_quadrature = geometry_S.volume();

      for (int k = 0; k < dimension; ++k) {
        innerValue += cutting_function * weight_micro_quadrature * diffusive_flux_y_S_inner[0][k] * unitOuterNormal[k];
        outerValue += cutting_function * weight_micro_quadrature * diffusive_flux_y_S_outer[0][k] * unitOuterNormal[k];
      }
    }

    local_indicator += pow(innerValue - outerValue, 2.0);

    // in 2D (and this is what we assume) it holds |E|*H_E^3 = H_E^4 = |E|^4
    return pow(edge_length, 4.0) * (1.0 / 2.0) * local_indicator;
  } // end of method

  // NOTE: This method ONLY works for uniformly refined micro-grid and GRIDDIM=2!
  // (even though it might be extended easily to other cases, here we use the cheapest strategy, which just works for
  // the case described above!)
  // here, uniform also means that we have the same number of elements in every direction!)
  // \bar{\eta}_T^{res}
  RangeType indicator_res_T(const EntityType& entity, const DiscreteFunctionType& u_H,
                            const PeriodicDiscreteFunctionType& corrector_u_H_on_entity) const {
    RangeType local_indicator(0.0);

    const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);

    const auto& globalEntityGeometry = entity.geometry();
    const auto& x_T = globalEntityGeometry.center();
    const auto& x_T_local = globalEntityGeometry.center();
    const auto entityVolume = globalEntityGeometry.volume();

    // \nabla u_H(x_T)
    auto u_H_local = u_H.localFunction(entity);
    JacobianRangeType gradient_u_H(0.);

    u_H_local.jacobian(x_T_local, gradient_u_H);

    if (dimension != 2) {
      std::cout << "The error indicator 'indicator_res_T' is not implemented for dimension!=2 and only works for "
                   "uniformly refined micro-grids!" << std::endl;
    }

    // edge length of a boundary face
    RangeType ref_edge_length = 1.0;
    const auto& auxGridPart = auxiliaryDiscreteFunctionSpace_.gridPart();
    // we just need one element to determine all the properties (due to uniform refinement)!
    for (const auto& intersection : DSC::intersectionRange(auxGridPart, *(auxiliaryDiscreteFunctionSpace_.begin()))) {
      const auto& faceGeometry = intersection.geometry();
      ref_edge_length = std::max(ref_edge_length, RangeType(faceGeometry.volume()));
    }

    // number of boundary faces per cube-edge:
    const int num_boundary_faces_per_direction = int((1 / ref_edge_length) + 0.2);
    // (+0.2 to avoid rounding errors)

    // generalized jump up/down
    std::vector<RangeType> jump_up_down(num_boundary_faces_per_direction);

    // generalized jump left/right
    std::vector<RangeType> jump_left_right(num_boundary_faces_per_direction);

    for (int id = 0; id < num_boundary_faces_per_direction; ++id) {
      jump_up_down[id] = 0.0;
      jump_left_right[id] = 0.0;
    }

    // procedure for computing the generalized jumps:
    // 1. compute the center of the current face E: (x_E,y_E)
    // 2. one and only one of these values fulfiles abs()=1/2
    // (i.e. we either have 'abs(x_E)=1/2' or 'abs(y_E)=1/2')
    // without loss of generality we assume 'abs(y_E)=1/2', then
    // we are in the setting of 'jump_up_down'
    // 3. Any 'jump_up_down' recieves a unique ID by
    // jump up down over (x,-1/2) and (x,1/2) = jump_up_down[ (num_boundary_faces_per_direction/2) + (( x -
    // (edge_length/2) ) / edge_length) ]
    // (num_boundary_faces_per_direction/2) + (( x - (edge_length/2) ) / edge_length) is the corresponding ID
    // 4. the situation is identical for 'abs(x_E)=1/2'.

    for (const auto& micro_entity : auxiliaryDiscreteFunctionSpace_) {
      const auto center_local = micro_entity.geometry().local(micro_entity.geometry().center());
      // Q_h(u_H)(x_T,y) on the micro entity:
      auto loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(micro_entity);
      JacobianRangeType gradient_Q_u_H_x_T(0.);
      loc_Q_u_H_x_T.jacobian(center_local, gradient_Q_u_H_x_T);

      // S denotes the micro grid element (i.e. 'micro_entity')
      const auto& geometry_S = micro_entity.geometry();
      // y_S denotes the barycenter of the micro grid element S:
      const auto& y_S = geometry_S.center();

      // to evaluate A^{\epsilon}_h (in center of current inner entity):
      DomainType globalPoint;
      JacobianRangeType direction_of_diffusion;
      for (int k = 0; k < dimension; ++k) {
        direction_of_diffusion[0][k] = gradient_u_H[0][k] + gradient_Q_u_H_x_T[0][k];
        globalPoint[k] = x_T[k] + (delta * y_S[k]);
      }

      JacobianRangeType diffusive_flux_y_S;
      diffusion_.diffusiveFlux(globalPoint, direction_of_diffusion, diffusive_flux_y_S);

      // ----------------------------------------------------

      for (const auto& intersection : DSC::intersectionRange(auxGridPart, entity)) {
        // Note: we are on the zero-centered unit cube! (That's why everything works!)
        const auto& faceGeometry = intersection.geometry();
        const auto edge_length = faceGeometry.volume();
        const auto unitOuterNormal = intersection.centerUnitOuterNormal();

        // if there is a neighbor entity (the normal gradient jumps)
        if (intersection.neighbor()) {
          // --------- the 'outer entity' (micro grid) ------------
          // ( neighbor entity )

          const auto outer_p_it = intersection.outside();
          const auto& outer_micro_entity = *outer_p_it;
          // S denotes the micro grid element (i.e. 'micro_entity')
          const auto& outer_geometry_S = outer_micro_entity.geometry();
          // outer_y_S denotes the barycenter of the neighbor micro grid element of S:
          const auto& outer_y_S = outer_geometry_S.center();
          const auto outer_y_S_local = outer_geometry_S.local(outer_y_S);

          // Q_h(u_H)(x_T,y) on the neighbor entity:
          auto outer_loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(outer_micro_entity);
          JacobianRangeType gradient_outer_Q_u_H_x_T(0.);
          outer_loc_Q_u_H_x_T.jacobian(outer_y_S_local, gradient_outer_Q_u_H_x_T);

          // to evaluate A^{\epsilon}_h (in center of current outer entity):
          DomainType outer_globalPoint;
          JacobianRangeType direction_of_diffusion_outside;
          for (int k = 0; k < dimension; ++k) {
            direction_of_diffusion_outside[0][k] = gradient_u_H[0][k] + gradient_outer_Q_u_H_x_T[0][k];
            outer_globalPoint[k] = x_T[k] + (delta * outer_y_S[k]);
          }

          JacobianRangeType diffusive_flux_y_S_outside;
          diffusion_.diffusiveFlux(outer_globalPoint, direction_of_diffusion_outside, diffusive_flux_y_S_outside);

          // ----------------------------------------------------

          RangeType aux_value = 0.0;

          for (int k = 0; k < dimension; ++k)
            aux_value += (diffusive_flux_y_S_outside[0][k] - diffusive_flux_y_S[0][k]) * unitOuterNormal[k];

          local_indicator += pow(edge_length, 4.0) * pow(aux_value, 2.0);
        } else {
          // if there is no neighbor entity, the face is a boundary face and we use the generalized gradient jumps.

          // Remember:
          // procedure for computing the generalized jumps:
          // 1. compute the center of the current face E: (x_E,y_E)
          // 2. one and only one of these values fulfiles abs()=1/2
          // (i.e. we either have 'abs(x_E)=1/2' or 'abs(y_E)=1/2')
          // without loss of generality we assume 'abs(y_E)=1/2', then
          // we are in the setting of 'jump_up_down'
          // 3. Any 'jump_up_down' recieves a unique ID by
          // jump up down over (x,-1/2) and (x,1/2) = ... jump_up_down[ (( x - (edge_length/2) ) / edge_length) ]
          // ... (( x - (edge_length/2) ) / edge_length) is the corresponding ID
          // 4. the situation is identical for 'abs(x_E)=1/2'.

          const auto& edge_center = geometry_S.center();

          if (fabs(edge_center[0]) == 0.5) {
            // + 0.2 to avoid rounding errors!
            const auto id = int((num_boundary_faces_per_direction / 2) +
                                ((edge_center[1] - (edge_length / 2.0)) / edge_length) + 0.2);

            // unit outer normal creates the correct sign!
            for (int k = 0; k < dimension; ++k)
              jump_up_down[id] += pow(edge_length, 2.0) * (diffusive_flux_y_S[0][k] * unitOuterNormal[k]);

            // std :: cout << "edge_center = " << edge_center << std :: endl;
            // std :: cout << "edge_length = " << edge_length << std :: endl;
            // std :: cout << "unitOuterNormal = " << unitOuterNormal << std :: endl;
            // std :: cout << "diffusive_flux_y_S = " << diffusive_flux_y_S[ 0 ] << std :: endl;
            // std :: cout << "jump_up_down id = " << id << std :: endl;
            // std :: cout << "jump_up_down[" << id << "] = " << jump_up_down[ id ] << std :: endl << std :: endl;
          } else {
            // + 0.2 to avoid rounding errors!
            int id = int((num_boundary_faces_per_direction / 2) +
                         ((edge_center[0] - (edge_length / 2.0)) / edge_length) + 0.2);

            // unit outer normal creates the correct sign!
            for (int k = 0; k < dimension; ++k)
              jump_left_right[id] += pow(edge_length, 2.0) * (diffusive_flux_y_S[0][k] * unitOuterNormal[k]);

            // std :: cout << "edge_center = " << edge_center << std :: endl;
            // std :: cout << "edge_length = " << edge_length << std :: endl;
            // std :: cout << "unitOuterNormal = " << unitOuterNormal << std :: endl;
            // std :: cout << "diffusive_flux_y_S = " << diffusive_flux_y_S[ 0 ] << std :: endl;
            // std :: cout << "jump_left_right id = " << id << std :: endl;
            // std :: cout << "jump_left_right[" << id << "] = " << jump_left_right[ id ] << std :: endl << std :: endl;
          }
        }
      }
    }

    for (int id = 0; id < num_boundary_faces_per_direction; ++id) {
      local_indicator += (pow(jump_up_down[id], 2.0) + pow(jump_left_right[id], 2.0));
    }

    local_indicator *= entityVolume;

    return local_indicator;
  } // end of method

  // only for TFR:

  // NOTE: This method ONLY works for uniformly refined micro-grid and GRIDDIM=2!
  // (even though it might be extended easily to other cases, here we use the cheapest strategy, which just works for
  // the case described above!)
  // here, uniform also means that we have the same number of elements in every direction!)
  // \eta_T^{tfr}
  RangeType indicator_tfr_1(const EntityType& entity, const DiscreteFunctionType& u_H,
                            const PeriodicDiscreteFunctionType& corrector_u_H_on_entity) const {
    RangeType local_indicator(0.0);

    const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);
    const double epsilon_estimated = DSC_CONFIG_GET("hmm.epsilon_guess", 1.0f);
    const auto& globalEntityGeometry = entity.geometry();
    const auto& x_T = globalEntityGeometry.center();
    const auto& x_T_local = globalEntityGeometry.local(x_T);
    const auto entityVolume = globalEntityGeometry.volume();
    // \nabla u_H(x_T)
    const auto u_H_local = u_H.localFunction(entity);
    JacobianRangeType gradient_u_H(0.);

    u_H_local.jacobian(x_T_local, gradient_u_H);

    if (dimension != 2) {
      std::cout << "The error indicator 'indicator_tfr_1' is not implemented for dimension!=2 and only works for "
                   "uniformly refined micro-grids!" << std::endl;
    }

    // edge length of a boundary face
    RangeType ref_edge_length = 1.0;
    const auto& auxGridPart = auxiliaryDiscreteFunctionSpace_.gridPart();
    // we just need one element to determine all the properties (due to uniform refinement)!
    for (const auto& intersection : DSC::intersectionRange(auxGridPart, entity)) {
      const auto& faceGeometry = intersection.geometry();
      ref_edge_length = std::max(ref_edge_length, RangeType(faceGeometry.volume()));
    }

    // number of boundary faces per (\epsilon/\delta-scaled) cube edge:
    const auto num_boundary_faces_per_direction = int(((epsilon_estimated / delta) / ref_edge_length) + 0.2);
    // (+0.2 to avoid rounding errors)

    // generalized jump up/down
    std::vector<RangeType> jump_up_down(num_boundary_faces_per_direction, RangeType(0.0));
    // generalized jump left/right
    std::vector<RangeType> jump_left_right(num_boundary_faces_per_direction, RangeType(0.0));

    // did you find a boundary edge of the \eps-\delta-cube? (you must find it!!)
    bool eps_delta_boundary_edge_found = false;

    // procedure for computing the generalized jumps:
    // 1. compute the center of the current face E: (x_E,y_E)
    // 2. one and only one of these values fulfiles abs()=1/2
    // (i.e. we either have 'abs(x_E)=1/2' or 'abs(y_E)=1/2')
    // without loss of generality we assume 'abs(y_E)=1/2', then
    // we are in the setting of 'jump_up_down'
    // 3. Any 'jump_up_down' recieves a unique ID by
    // jump up down over (x,-1/2) and (x,1/2) = jump_up_down[ (num_boundary_faces_per_direction/2) + (( x -
    // (edge_length/2) ) / edge_length) ]
    // (num_boundary_faces_per_direction/2) + (( x - (edge_length/2) ) / edge_length) is the corresponding ID
    // 4. the situation is identical for 'abs(x_E)=1/2'.

    // iterator over the elements of the periodic micro grid:
    for (const auto& micro_entity : auxiliaryDiscreteFunctionSpace_) {
      // S denotes the micro grid element (i.e. 'micro_entity')
      const auto& geometry_S = micro_entity.geometry();
      // y_S denotes the barycenter of the micro grid element S:
      const auto& y_S = geometry_S.center();
      const auto& y_S_local = geometry_S.local(y_S);

      // Q_h(u_H)(x_T,y) on the micro entity:
      auto loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(micro_entity);
      JacobianRangeType gradient_Q_u_H_x_T(0.);
      loc_Q_u_H_x_T.jacobian(y_S_local, gradient_Q_u_H_x_T);

      // to evaluate A^{\epsilon}_h (in center of current inner entity):
      DomainType globalPoint;
      JacobianRangeType direction_of_diffusion;
      for (int k = 0; k < dimension; ++k) {
        direction_of_diffusion[0][k] = gradient_u_H[0][k] + gradient_Q_u_H_x_T[0][k];
        globalPoint[k] = (delta * y_S[k]) + x_T[k];
      }

      JacobianRangeType diffusive_flux_y_S;
      diffusion_.diffusiveFlux(globalPoint, direction_of_diffusion, diffusive_flux_y_S);

      if ((fabs(y_S[1]) <= ((0.5 * (epsilon_estimated / delta)))) &&
          (fabs(y_S[0]) <= ((0.5 * (epsilon_estimated / delta))))) {

        for (const auto& intersection : DSC::intersectionRange(auxGridPart, entity)) {
          // Note: we are on the zero-centered unit cube! (That's why everything works!)

          const auto& faceGeometry = intersection.geometry();
          const auto edge_length = faceGeometry.volume();
          const auto unitOuterNormal = intersection.centerUnitOuterNormal();
          const auto& edge_center = geometry_S.center();

          if ((fabs(edge_center[0]) == ((0.5 * (epsilon_estimated / delta)))) ||
              (fabs(edge_center[1]) == ((0.5 * (epsilon_estimated / delta))))) {
            // we use the generalized gradient jumps.

            if ((fabs(edge_center[0]) == ((0.5 * (epsilon_estimated / delta)))) &&
                (fabs(edge_center[1]) < ((0.5 * (epsilon_estimated / delta))))) {
              // + 0.2 to avoid rounding errors!
              int id = int((num_boundary_faces_per_direction / 2) +
                           ((edge_center[1] - (edge_length / 2.0)) / edge_length) + 0.2);

              // unit outer normal creates the correct sign!
              for (int k = 0; k < dimension; ++k)
                jump_up_down[id] += sqrt(edge_length) * (diffusive_flux_y_S[0][k] * unitOuterNormal[k]);

              // std :: cout << "edge_center = " << edge_center << std :: endl;
              // std :: cout << "edge_length = " << edge_length << std :: endl;
              // std :: cout << "unitOuterNormal = " << unitOuterNormal << std :: endl;
              // std :: cout << "diffusive_flux_y_S = " << diffusive_flux_y_S[ 0 ] << std :: endl;
              // std :: cout << "jump_up_down id = " << id << std :: endl;
              // std :: cout << "jump_up_down[" << id << "] = " << jump_up_down[ id ] << std :: endl << std :: endl;

              eps_delta_boundary_edge_found = true;
              // std :: cout << "Found eps/delta boundary edge with center = " << edge_center << std :: endl;
            }

            if ((fabs(edge_center[1]) == ((0.5 * (epsilon_estimated / delta)))) &&
                (fabs(edge_center[0]) < ((0.5 * (epsilon_estimated / delta))))) {
              // + 0.2 to avoid rounding errors!
              int id = int((num_boundary_faces_per_direction / 2) +
                           ((edge_center[0] - (edge_length / 2.0)) / edge_length) + 0.2);

              // unit outer normal creates the correct sign!
              for (int k = 0; k < dimension; ++k)
                jump_left_right[id] += sqrt(edge_length) * (diffusive_flux_y_S[0][k] * unitOuterNormal[k]);

              // std :: cout << "edge_center = " << edge_center << std :: endl;
              // std :: cout << "edge_length = " << edge_length << std :: endl;
              // std :: cout << "unitOuterNormal = " << unitOuterNormal << std :: endl;
              // std :: cout << "diffusive_flux_y_S = " << diffusive_flux_y_S[ 0 ] << std :: endl;
              // std :: cout << "jump_left_right id = " << id << std :: endl;
              // std :: cout << "jump_left_right[" << id << "] = " << jump_left_right[ id ] << std :: endl << std ::
              // endl;

              eps_delta_boundary_edge_found = true;
              // std :: cout << "Found eps/delta boundary edge with center = " << edge_center << std :: endl;
            }
          }

          if ((fabs(edge_center[0]) < (0.5 * (epsilon_estimated / delta))) &&
              (fabs(edge_center[1]) < (0.5 * (epsilon_estimated / delta)))) {
            // in this situation, p_nit always has an outside neighbor entity
            // (use the normal gradient jumps)

            // --------- the 'outer entity' (micro grid) ------------
            // ( neighbor entity )

            const auto outer_p_it = intersection.outside();
            const auto& outer_micro_entity = *outer_p_it;
            // S denotes the micro grid element (i.e. 'micro_entity')
            const auto& outer_geometry_S = outer_micro_entity.geometry();
            // outer_y_S denotes the barycenter of the neighbor micro grid element of S:
            const auto& outer_y_S = outer_geometry_S.center();
            const auto& outer_y_S_local = outer_geometry_S.local(outer_y_S);

            // Q_h(u_H)(x_T,y) on the neighbor entity:
            const auto outer_loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(outer_micro_entity);
            JacobianRangeType gradient_outer_Q_u_H_x_T(0.);
            outer_loc_Q_u_H_x_T.jacobian(outer_y_S_local, gradient_outer_Q_u_H_x_T);

            // to evaluate A^{\epsilon}_h (in center of current outer entity):
            DomainType outer_globalPoint;
            JacobianRangeType direction_of_diffusion_outside;
            for (int k = 0; k < dimension; ++k) {
              outer_globalPoint[k] = x_T[k] + (delta * outer_y_S[k]);
              direction_of_diffusion_outside[0][k] = gradient_u_H[0][k] + gradient_outer_Q_u_H_x_T[0][k];
            }

            JacobianRangeType diffusive_flux_y_S_outside;
            diffusion_.diffusiveFlux(outer_globalPoint, direction_of_diffusion_outside, diffusive_flux_y_S_outside);

            // ----------------------------------------------------

            RangeType aux_value = 0.0;

            for (int k = 0; k < dimension; ++k)
              aux_value += (diffusive_flux_y_S_outside[0][k] - diffusive_flux_y_S[0][k]) * unitOuterNormal[k];

            local_indicator += edge_length * pow(aux_value, 2.0);
          }
        } // end Intersection iterator
      }   // endif: "y_S \in \eps/\delta Y"
    }

    for (int id = 0; id < num_boundary_faces_per_direction; ++id) {
      local_indicator += (pow(jump_up_down[id], 2.0) + pow(jump_left_right[id], 2.0));
    }

    if (!eps_delta_boundary_edge_found) {
      DUNE_THROW(Dune::InvalidStateException, "Error! Make sure that the restriction of the Y-triangulation on "
                                              "'eps/delta Y' is a complete periodic triangulation of 'eps/delta Y' on "
                                              "its own. (for instance: delta = 2 epsilon should work)");
    }

    local_indicator *= entityVolume;

    return local_indicator;
  } // end of method

  // NOTE: This method ONLY works for uniformly refined micro-grid and GRIDDIM=2!
  // (even though it might be extended easily to other cases, here we use the cheapest strategy, which just works for
  // the case described above!)
  // here, uniform also means that we have the same number of elements in every direction!)
  // This indicator does not exist in theory. It is is additionally multiplied with h^3 to fit the other orders of
  // convergence. In this setting it is more suitable for a comparison to capture the effect of boundary jumps in the
  // case of a wrong boundary condition
  RangeType indicator_effective_tfr(const EntityType& entity, const DiscreteFunctionType& u_H,
                                    const PeriodicDiscreteFunctionType& corrector_u_H_on_entity) const {
    RangeType local_indicator(0.0);

    const double delta = DSC_CONFIG_GET("hmm.delta", 1.0f);
    const double epsilon_estimated = DSC_CONFIG_GET("hmm.epsilon_guess", 1.0f);

    const auto& globalEntityGeometry = entity.geometry();
    const auto& x_T = globalEntityGeometry.center();
    const auto& x_T_local = globalEntityGeometry.local(x_T);
    const auto entityVolume = globalEntityGeometry.volume();

    // \nabla u_H(x_T)
    auto u_H_local = u_H.localFunction(entity);
    JacobianRangeType gradient_u_H(0.);

    u_H_local.jacobian(x_T_local, gradient_u_H);

    static_assert(dimension != 2, "The error indicator 'indicator_tfr_1' is not implemented for dimension!=2 and only "
                                  "works for uniformly refined micro-grids!");

    // edge length of a boundary face
    RangeType ref_edge_length(0.);
    const auto& auxGridPart = auxiliaryDiscreteFunctionSpace_.gridPart();
    // we just need one element to determine all the properties (due to uniform refinement)!
    for (const auto& intersection : DSC::intersectionRange(auxGridPart, entity)) {
      const auto& faceGeometry = intersection.geometry();
      ref_edge_length = std::max(ref_edge_length, RangeType(faceGeometry.volume()));
    }

    // number of boundary faces per (\epsilon/\delta-scaled) cube edge:
    const int num_boundary_faces_per_direction = int(((epsilon_estimated / delta) / ref_edge_length) + 0.2);
    // (+0.2 to avoid rounding errors)

    // generalized jump up/down
    std::vector<RangeType> jump_up_down(num_boundary_faces_per_direction, RangeType(0.0));
    // generalized jump left/right
    std::vector<RangeType> jump_left_right(num_boundary_faces_per_direction, RangeType(0.0));

    // did you find a boundary edge of the \eps-\delta-cube? (you must find it!!)
    bool eps_delta_boundary_edge_found = false;

    // procedure for computing the generalized jumps:
    // 1. compute the center of the current face E: (x_E,y_E)
    // 2. one and only one of these values fulfiles abs()=1/2
    // (i.e. we either have 'abs(x_E)=1/2' or 'abs(y_E)=1/2')
    // without loss of generality we assume 'abs(y_E)=1/2', then
    // we are in the setting of 'jump_up_down'
    // 3. Any 'jump_up_down' recieves a unique ID by
    // jump up down over (x,-1/2) and (x,1/2) = jump_up_down[ (num_boundary_faces_per_direction/2) + (( x -
    // (edge_length/2) ) / edge_length) ]
    // (num_boundary_faces_per_direction/2) + (( x - (edge_length/2) ) / edge_length) is the corresponding ID
    // 4. the situation is identical for 'abs(x_E)=1/2'.

    for (const auto& micro_entity : auxiliaryDiscreteFunctionSpace_) {
      // S denotes the micro grid element (i.e. 'micro_entity')
      const auto& geometry_S = micro_entity.geometry();
      // y_S denotes the barycenter of the micro grid element S:
      const auto& y_S = geometry_S.center();
      const auto& y_S_local = geometry_S.local(y_S);
      // Q_h(u_H)(x_T,y) on the micro entity:
      auto loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(micro_entity);
      JacobianRangeType gradient_Q_u_H_x_T(0.);
      loc_Q_u_H_x_T.jacobian(y_S_local, gradient_Q_u_H_x_T);
      // to evaluate A^{\epsilon}_h (in center of current inner entity):
      DomainType globalPoint;
      JacobianRangeType direction_of_diffusion;
      for (int k = 0; k < dimension; ++k) {
        direction_of_diffusion[0][k] = gradient_u_H[0][k] + gradient_Q_u_H_x_T[0][k];
        globalPoint[k] = x_T[k] + (delta * y_S[k]);
      }

      JacobianRangeType diffusive_flux_y_S;
      diffusion_.diffusiveFlux(globalPoint, direction_of_diffusion, diffusive_flux_y_S);

      if ((fabs(y_S[1]) <= ((0.5 * (epsilon_estimated / delta)))) &&
          (fabs(y_S[0]) <= ((0.5 * (epsilon_estimated / delta))))) {

        for (const auto& intersection : DSC::intersectionRange(auxGridPart, micro_entity)) {
          // Note: we are on the zero-centered unit cube! (That's why everything works!)

          const auto& faceGeometry = intersection.geometry();
          const auto edge_length = faceGeometry.volume();
          const auto unitOuterNormal = intersection.centerUnitOuterNormal();
          const auto& edge_center = geometry_S.center();

          if ((fabs(edge_center[0]) == ((0.5 * (epsilon_estimated / delta)))) ||
              (fabs(edge_center[1]) == ((0.5 * (epsilon_estimated / delta))))) {
            // we use the generalized gradient jumps.

            if ((fabs(edge_center[0]) == ((0.5 * (epsilon_estimated / delta)))) &&
                (fabs(edge_center[1]) < ((0.5 * (epsilon_estimated / delta))))) {
              // + 0.2 to avoid rounding errors!
              int id = int((num_boundary_faces_per_direction / 2) +
                           ((edge_center[1] - (edge_length / 2.0)) / edge_length) + 0.2);

              // unit outer normal creates the correct sign!
              for (int k = 0; k < dimension; ++k)
                jump_up_down[id] += pow(edge_length, 2.0) * (diffusive_flux_y_S[0][k] * unitOuterNormal[k]);

              // std :: cout << "edge_center = " << edge_center << std :: endl;
              // std :: cout << "edge_length = " << edge_length << std :: endl;
              // std :: cout << "unitOuterNormal = " << unitOuterNormal << std :: endl;
              // std :: cout << "diffusive_flux_y_S = " << diffusive_flux_y_S[ 0 ] << std :: endl;
              // std :: cout << "jump_up_down id = " << id << std :: endl;
              // std :: cout << "jump_up_down[" << id << "] = " << jump_up_down[ id ] << std :: endl << std :: endl;

              eps_delta_boundary_edge_found = true;
              // std :: cout << "Found eps/delta boundary edge with center = " << edge_center << std :: endl;
            }

            if ((fabs(edge_center[1]) == ((0.5 * (epsilon_estimated / delta)))) &&
                (fabs(edge_center[0]) < ((0.5 * (epsilon_estimated / delta))))) {
              // + 0.2 to avoid rounding errors!
              int id = int((num_boundary_faces_per_direction / 2) +
                           ((edge_center[0] - (edge_length / 2.0)) / edge_length) + 0.2);

              // unit outer normal creates the correct sign!
              for (int k = 0; k < dimension; ++k)
                jump_left_right[id] += pow(edge_length, 2.0) * (diffusive_flux_y_S[0][k] * unitOuterNormal[k]);

              // std :: cout << "edge_center = " << edge_center << std :: endl;
              // std :: cout << "edge_length = " << edge_length << std :: endl;
              // std :: cout << "unitOuterNormal = " << unitOuterNormal << std :: endl;
              // std :: cout << "diffusive_flux_y_S = " << diffusive_flux_y_S[ 0 ] << std :: endl;
              // std :: cout << "jump_left_right id = " << id << std :: endl;
              // std :: cout << "jump_left_right[" << id << "] = " << jump_left_right[ id ] << std :: endl << std ::
              // endl;

              eps_delta_boundary_edge_found = true;
              // std :: cout << "Found eps/delta boundary edge with center = " << edge_center << std :: endl;
            }
          }

          if ((fabs(edge_center[0]) < (0.5 * (epsilon_estimated / delta))) &&
              (fabs(edge_center[1]) < (0.5 * (epsilon_estimated / delta)))) {
            // in this situation, p_nit always has an outside neighbor entity
            // (use the normal gradient jumps)

            // --------- the 'outer entity' (micro grid) ------------
            // ( neighbor entity )

            const auto outer_p_it = intersection.outside();
            const auto& outer_micro_entity = *outer_p_it;
            // S denotes the micro grid element (i.e. 'micro_entity')
            const auto& outer_geometry_S = outer_micro_entity.geometry();
            // outer_y_S denotes the barycenter of the neighbor micro grid element of S:
            const auto& outer_y_S = outer_geometry_S.center();
            const auto& outer_y_S_local = outer_geometry_S.local(outer_y_S);

            // Q_h(u_H)(x_T,y) on the neighbor entity:
            auto outer_loc_Q_u_H_x_T = corrector_u_H_on_entity.localFunction(outer_micro_entity);
            JacobianRangeType gradient_outer_Q_u_H_x_T(0.);
            outer_loc_Q_u_H_x_T.jacobian(outer_y_S_local, gradient_outer_Q_u_H_x_T);

            // to evaluate A^{\epsilon}_h (in center of current outer entity):
            DomainType outer_globalPoint;
            JacobianRangeType direction_of_diffusion_outside;
            for (int k = 0; k < dimension; ++k) {
              outer_globalPoint[k] = x_T[k] + (delta * outer_y_S[k]);
              direction_of_diffusion_outside[0][k] = gradient_u_H[0][k] + gradient_outer_Q_u_H_x_T[0][k];
            }

            JacobianRangeType diffusive_flux_y_S_outside;
            diffusion_.diffusiveFlux(outer_globalPoint, direction_of_diffusion_outside, diffusive_flux_y_S_outside);

            // ----------------------------------------------------

            RangeType aux_value = 0.0;

            for (int k = 0; k < dimension; ++k)
              aux_value += (diffusive_flux_y_S_outside[0][k] - diffusive_flux_y_S[0][k]) * unitOuterNormal[k];

            local_indicator += pow(edge_length, 4.0) * pow(aux_value, 2.0);
          }
        } // end Intersection iterator
      }   // endif: "y_S \in \eps/\delta Y"
    }

    for (int id = 0; id < num_boundary_faces_per_direction; ++id) {
      local_indicator += (pow(jump_up_down[id], 2.0) + pow(jump_left_right[id], 2.0));
    }

    if (!eps_delta_boundary_edge_found) {
      DUNE_THROW(Dune::InvalidStateException, "Error! Make sure that the restriction of the Y-triangulation on "
                                              "'eps/delta Y' is a complete periodic triangulation of 'eps/delta Y' on "
                                              "its own. (for instance: delta = 2 epsilon should work)");
    }

    local_indicator *= entityVolume;

    return local_indicator;
  } // end of method
};  // end of class ErrorEstimater

} // namespace HMM {
} // namespace Multiscale {
} // namespace Dune {

#endif // ifndef DUNE_ERRORESTIMATER_HH
