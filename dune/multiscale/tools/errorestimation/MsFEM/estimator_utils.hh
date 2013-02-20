#ifndef ESTIMATOR_UTILS_HH
#define ESTIMATOR_UTILS_HH

#include <dune/multiscale/tools/misc.hh>

namespace Dune {
namespace Multiscale {
namespace MsFEM {

template < class EstimatorType >
struct EstimatorUtils {
  //! create N hostgrid functions from N subgridfunctions
  template < std::array<int, 1>::size_type N >
  static void subgrid_to_hostrid_function(const std::array<typename EstimatorType::SubGridDiscreteFunctionType, N>& sub_funcs,
                                          std::array< typename EstimatorType::DiscreteFunctionPointer, N >& host_funcs) {
    for( auto& host_func: host_funcs)
      host_func->clear();

    const auto& subGrid = sub_funcs[0].space().grid();
    for (const auto& sub_entity : sub_funcs[0].space())
    {
      const auto host_entity_pointer = subGrid.template getHostEntity< 0 >(sub_entity);
      const auto& host_entity = *host_entity_pointer;
      for (std::array<int, 1>::size_type j = 0; j < N; ++j) {
        const auto sub_loc_value = sub_funcs[j].localFunction(sub_entity);
        auto host_loc_value = host_funcs[j]->localFunction(host_entity);
        const auto numBaseFunctions = sub_loc_value.baseFunctionSet().size();
        for (unsigned int i = 0; i < numBaseFunctions; ++i)
        {
          host_loc_value[i] = sub_loc_value[i];
        }
      }
    }
  } // subgrid_to_hostrid_function


  //! is a given point on a given face?
  static bool point_on_face(const typename EstimatorType::Intersection& face,
                            const typename EstimatorType::DomainType& point) {
    const auto corner_0 = face.geometry().corner(0);
    const auto corner_1 = face.geometry().corner(1);

    if (corner_0[0] == corner_1[0])
    {
      if (point[0] != corner_0[0]) {
        return false;
      } else {
        //!TODO das face ist nen punkt ungleich point, warum sollte rückgabe true noch möglich sein?
        const typename EstimatorType::RangeType lambda = (point[1] - corner_1[1])
                                                         / (corner_0[1] - corner_1[1]);
        return ( (lambda >= 0.0) && (lambda <= 1.0) );
      }
    } else {
      const typename EstimatorType::RangeType lambda = (point[0] - corner_1[0])
                                                       / (corner_0[0] - corner_1[0]);

      if ( (lambda >= 0.0) && (lambda <= 1.0) )
      {
        const typename EstimatorType::RangeType convex_comb = (lambda * corner_0[1])
                                                              + ( (1.0 - lambda) * corner_1[1] );
        return convex_comb == point[1];
      } else {
        return false;
      }
    }
  } // point_on_face

  // is a given face part of another given coarse face
  static bool is_subface(const typename EstimatorType::Intersection& fine_face,
                         const typename EstimatorType::Intersection& coarse_face) {
    const auto corner_0 = fine_face.geometry().corner(0);
    const auto corner_1 = fine_face.geometry().corner(1);
    return ( point_on_face(coarse_face, corner_0) && point_on_face(coarse_face, corner_1) );
  } // is_subface


  static std::pair<typename EstimatorType::JumpArray,typename EstimatorType::JumpArray>
  flux_contributions( const typename EstimatorType::SubGridDiscreteFunctionSpaceType& localDiscreteFunctionSpace,
                      const typename EstimatorType::SubGridType& sub_grid_U_T,
                      const typename EstimatorType::LeafIndexSetType& coarseGridLeafIndexSet,
                      const typename EstimatorType::DiscreteFunctionPointerPair& cflux_coarse_ent_host,
                      const typename EstimatorType::DiscreteFunctionType& msfem_coarse_part,
                      typename EstimatorType::FluxContainerType& cflux_neighbor_ent_host,
                      const int index_coarse_entity,
                      const std::array<typename EstimatorType::RangeType, 3>& coarse_face_volume,
                      const int level_difference,
                      const typename EstimatorType::DiscreteFunctionSpaceType& fineDiscreteFunctionSpace)
  {
    // jump for each face
    typename EstimatorType::JumpArray jump = {{ 0.0, 0.0, 0.0 }};
    // coarse grid jump for each face
    typename EstimatorType::JumpArray coarse_jump = {{ 0.0, 0.0, 0.0 }};

    for (const auto& sub_entity : localDiscreteFunctionSpace)
    {
      auto host_entity_pointer = sub_grid_U_T.template getHostEntity< 0 >(sub_entity);
      const auto& host_entity = *host_entity_pointer;

      auto father_of_sub_grid_entity = Stuff::Grid::make_father(coarseGridLeafIndexSet,
                                                                host_entity_pointer,
                                                                level_difference);
      const int coarse_sub_father_index = coarseGridLeafIndexSet.index(*father_of_sub_grid_entity);
      if (coarse_sub_father_index != index_coarse_entity)
      { continue; }

      const auto loc_cf_coarse_ent_e0 = cflux_coarse_ent_host[0]->localFunction(host_entity);
      const auto loc_cf_coarse_ent_e1 = cflux_coarse_ent_host[1]->localFunction(host_entity);

      const auto loc_msfem_coarse_part = msfem_coarse_part.localFunction(host_entity);

      auto end_it_U_T = fineDiscreteFunctionSpace.gridPart().iend(host_entity);
      for (auto face_it_U_T = fineDiscreteFunctionSpace.gridPart().ibegin(host_entity);
           face_it_U_T != end_it_U_T; ++face_it_U_T)
      {
        int relevant_face_index = cflux_neighbor_ent_host.intersection_compatible(*face_it_U_T);

        if ( (relevant_face_index == -1) || (!face_it_U_T->neighbor()) )
        { continue; }

        auto outside_sub_it = face_it_U_T->outside();
        assert(cflux_neighbor_ent_host.fluxes[relevant_face_index][0]);
        assert(cflux_neighbor_ent_host.fluxes[relevant_face_index][1]);
        auto loc_cf_coarse_neighbor_ent_e0
          = (*cflux_neighbor_ent_host.fluxes[relevant_face_index][0]).localFunction(host_entity);
        auto loc_cf_coarse_neighbor_ent_e1
          = (*cflux_neighbor_ent_host.fluxes[relevant_face_index][1]).localFunction(host_entity);

        auto loc_msfem_coarse_part_neighbor = msfem_coarse_part.localFunction(*outside_sub_it);

        // evaluate the gradient of the MsfEM coarse part in the center of the coarse entity
        const typename EstimatorType::EntityQuadratureType coarseEntQuadrature(host_entity, 0);
        typename EstimatorType::JacobianRangeType gradient_msfem_coarse_ent(0.);
        loc_msfem_coarse_part.jacobian(coarseEntQuadrature[0], gradient_msfem_coarse_ent);

        // evaluate the gradient of the MsfEM coarse part in the center of the current neighbor of the coarse entity
        const typename EstimatorType::EntityQuadratureType coarseEntQuadratureNeighbor(*outside_sub_it, 0);
        typename EstimatorType::JacobianRangeType gradient_msfem_coarse_neighbor_ent(0.);
        loc_msfem_coarse_part_neighbor.jacobian(coarseEntQuadratureNeighbor[0], gradient_msfem_coarse_neighbor_ent);

        const typename EstimatorType::FaceQuadratureType faceQuadrature(
          fineDiscreteFunctionSpace.gridPart(), *face_it_U_T,
          2 * fineDiscreteFunctionSpace.order() + 2, EstimatorType::FaceQuadratureType::INSIDE);
        // inside macht hier keinen Unterschied, da wir formal stetige Funktionen haben und nicht die Gradienten
        // auswerten

        const auto& faceGeometry = face_it_U_T->geometry();
        const std::size_t numQuadraturePoints = faceQuadrature.nop();

        typedef typename EstimatorType::RangeType RangeType;
        RangeType check_sum(0.0);
        for (std::size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
        {
          const auto local_point = faceGeometry.local( faceQuadrature.point(quadraturePoint) );

          // integration factors
          const double integrationFactor = faceGeometry.integrationElement(local_point);
          // weight
          const double quadratureWeight = faceQuadrature.weight(quadraturePoint);

          check_sum += integrationFactor * quadratureWeight;

          RangeType value_ent_e0;
          loc_cf_coarse_ent_e0.evaluate(faceQuadrature[quadraturePoint], value_ent_e0);

          RangeType value_neighbor_ent_e0;
          loc_cf_coarse_neighbor_ent_e0.evaluate(faceQuadrature[quadraturePoint], value_neighbor_ent_e0);

          RangeType value_ent_e1;
          loc_cf_coarse_ent_e1.evaluate(faceQuadrature[quadraturePoint], value_ent_e1);

          RangeType value_neighbor_ent_e1;
          loc_cf_coarse_neighbor_ent_e1.evaluate(faceQuadrature[quadraturePoint], value_neighbor_ent_e1);

          RangeType H_E = coarse_face_volume[relevant_face_index];

          // wenn man tunen moechtet... fabs() einbinden und das Vorzeichen davor aendern.

          RangeType jump_contribution = gradient_msfem_coarse_ent[0][0]
                                        * ( /*fabs*/ (value_ent_e0) + /*-fabs*/ (value_neighbor_ent_e0) );
          jump_contribution += gradient_msfem_coarse_ent[0][1]
                               * ( /*fabs*/ (value_ent_e1) + /*-fabs*/ (value_neighbor_ent_e1) );

          jump[relevant_face_index] += H_E * integrationFactor * quadratureWeight * pow(jump_contribution, 2.0);

          jump_contribution = gradient_msfem_coarse_neighbor_ent[0][0]
                              * ( /*fabs*/ (value_ent_e0) + /*-fabs*/ (value_neighbor_ent_e0) );
          jump_contribution += gradient_msfem_coarse_neighbor_ent[0][1]
                               * ( /*fabs*/ (value_ent_e1) + /*-fabs*/ (value_neighbor_ent_e1) );

          jump[relevant_face_index] += H_E * integrationFactor * quadratureWeight * pow(jump_contribution, 2.0);

          typename EstimatorType::JacobianRangeType coarse_jump_contribution;
          coarse_jump_contribution[0] = gradient_msfem_coarse_ent[0] - gradient_msfem_coarse_neighbor_ent[0];

          RangeType cjump(0.0);
          cjump += fabs(value_ent_e0 * coarse_jump_contribution[0][0] + value_ent_e1 * coarse_jump_contribution[0][1]);
          cjump += fabs(
            value_neighbor_ent_e0 * coarse_jump_contribution[0][0] + value_neighbor_ent_e1
            * coarse_jump_contribution[0][1]);

          coarse_jump[relevant_face_index] += H_E * integrationFactor * quadratureWeight * pow(cjump, 2.0);

        } // done loop over all quadrature points

        if ( check_sum != faceGeometry.volume() )
        {
          DUNE_THROW(Dune::InvalidStateException, "Error in Face Quadrature.");
        }
      }
    }
    return std::make_pair(jump, coarse_jump);
  }

};

} //namespace MsFEM {
} //namespace Multiscale {
} //namespace Dune {

#endif // ESTIMATOR_UTILS_HH
