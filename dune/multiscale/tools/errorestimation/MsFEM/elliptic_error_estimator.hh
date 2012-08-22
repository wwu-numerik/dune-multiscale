#ifndef DUNE_MSFEM_ERRORESTIMATOR_HH
#define DUNE_MSFEM_ERRORESTIMATOR_HH

#include <dune/common/unused.hh>

// where the quadratures are defined
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/multiscale/tools/errorestimation/MsFEM/conservative_flux_solver.hh>

namespace Dune {
template< class DiscreteFunctionImp,
          class DiffusionImp,
          class SourceImp,
          class MacroMicroGridSpecifierImp,
          class SubGridListImp >
class MsFEMErrorEstimator
{
  typedef DiffusionImp               DiffusionOperatorType;
  typedef SourceImp                  SourceType;
  typedef MacroMicroGridSpecifierImp MacroMicroGridSpecifierType;
  typedef SubGridListImp             SubGridListType;

  //! Necessary typedefs for the DiscreteFunctionImp:

  typedef DiscreteFunctionImp                              DiscreteFunctionType;
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;

  typedef typename DiscreteFunctionType::LocalFunctionType         LocalFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::RangeFieldType    RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::DomainType        DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType         RangeType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionType::DofIteratorType           DofIteratorType;
  typedef typename DiscreteFunctionSpaceType::GridPartType         GridPartType;
  typedef typename DiscreteFunctionSpaceType::GridType             GridType;
  typedef typename DiscreteFunctionSpaceType::IteratorType         IteratorType;
  typedef typename DiscreteFunctionSpaceType::LagrangePointSetType LagrangePointSetType;

  typedef typename GridType::Traits::LeafIndexSet LeafIndexSetType;

  typedef typename GridPartType::IntersectionIteratorType         IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection         Intersection;
  typedef typename GridType::template Codim< 0 >::Entity          EntityType;
  typedef typename GridType::template Codim< 0 >::EntityPointer   EntityPointerType;
  typedef typename GridType::template Codim< 1 >::Entity          FaceType;
  typedef typename GridType::template Codim< 1 >::EntityPointer   FacePointerType;
  typedef typename GridType::template Codim< 0 >::Geometry        EntityGeometryType;
  typedef typename GridType::template Codim< 1 >::Geometry        FaceGeometryType;
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

  typedef CachingQuadrature< GridPartType, 0 > EntityQuadratureType;
  typedef CachingQuadrature< GridPartType, 1 > FaceQuadratureType;

  // --------------------------- subgrid typedefs ------------------------------------

  typedef SubGrid< GridType::dimension, GridType > SubGridType;
  typedef LeafGridPart< SubGridType >              SubGridPart;

  typedef typename SubGridType::Traits::LeafIndexSet SubGridLeafIndexSet;

  typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, SubGridPart,
                                         1 /*POLORDER*/ > SubGridDiscreteFunctionSpaceType;

  typedef AdaptiveDiscreteFunction< SubGridDiscreteFunctionSpaceType > SubGridDiscreteFunctionType;

  typedef typename SubGridDiscreteFunctionSpaceType::IteratorType SubGridIteratorType;

  typedef typename SubGridIteratorType::Entity SubGridEntityType;

  typedef typename SubGridEntityType::EntityPointer SubGridEntityPointerType;

  typedef typename SubGridDiscreteFunctionType::LocalFunctionType SubGridLocalFunctionType;

  typedef typename SubGridDiscreteFunctionSpaceType::LagrangePointSetType
  SubGridLagrangePointSetType;

  enum { faceCodim = 1 };
  typedef typename SubGridLagrangePointSetType::template Codim< faceCodim >
    ::SubEntityIteratorType
  SubGridFaceDofIteratorType;

  typedef typename LagrangePointSetType::template Codim< faceCodim >
    ::SubEntityIteratorType
  FaceDofIteratorType;

  typedef typename SubGridType::template Codim< 0 >::Geometry SubGridEntityGeometryType;

  typedef CachingQuadrature< SubGridPart, 0 > LocalGridEntityQuadratureType;

  //!-----------------------------------------------------------------------------------------

  enum { dimension = GridType::dimension };
  enum { spacePolOrd = DiscreteFunctionSpaceType::polynomialOrder };
  enum { maxnumOfBaseFct = 100 };

  const DiscreteFunctionSpaceType& fineDiscreteFunctionSpace_;
  MacroMicroGridSpecifierType& specifier_;
  SubGridListType& subgrid_list_;
  const DiffusionOperatorType& diffusion_;
  const SourceType& f_;
  std::string& path_;

public:
  MsFEMErrorEstimator(const DiscreteFunctionSpaceType& fineDiscreteFunctionSpace,
                      MacroMicroGridSpecifierType& specifier,
                      SubGridListType& subgrid_list,
                      const DiffusionOperatorType& diffusion,
                      const SourceType& f,
                      std::string& path)
    : fineDiscreteFunctionSpace_(fineDiscreteFunctionSpace)
      , specifier_(specifier)
      , subgrid_list_(subgrid_list)
      , diffusion_(diffusion)
      , f_(f)
      , path_(path)
  {}

  // create a hostgrid function from a subgridfunction
  void subgrid_to_hostrid_function(const SubGridDiscreteFunctionType& sub_func,
                                   DiscreteFunctionType& host_func) {
    host_func.clear();

    const SubGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace = sub_func.space();
    const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

    SubGridIteratorType sub_endit = subDiscreteFunctionSpace.end();
    for (SubGridIteratorType sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
    {
      const SubGridEntityType& sub_entity = *sub_it;

      EntityPointerType host_entity_pointer = subGrid.template getHostEntity< 0 >(*sub_it);
      const EntityType& host_entity = *host_entity_pointer;

      SubGridLocalFunctionType sub_loc_value = sub_func.localFunction(sub_entity);
      LocalFunctionType host_loc_value = host_func.localFunction(host_entity);

      const unsigned int numBaseFunctions = sub_loc_value.baseFunctionSet().numBaseFunctions();
      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        host_loc_value[i] = sub_loc_value[i];
      }
    }
  } // subgrid_to_hostrid_function

  // create twoa hostgrid functions from two subgridfunctions
  void subgrid_to_hostrid_function(const SubGridDiscreteFunctionType& sub_func_1,
                                   const SubGridDiscreteFunctionType& sub_func_2,
                                   DiscreteFunctionType& host_func_1,
                                   DiscreteFunctionType& host_func_2) {
    host_func_1.clear();
    host_func_2.clear();

    const SubGridDiscreteFunctionSpaceType& subDiscreteFunctionSpace = sub_func_1.space();
    const SubGridType& subGrid = subDiscreteFunctionSpace.grid();

    SubGridIteratorType sub_endit = subDiscreteFunctionSpace.end();
    for (SubGridIteratorType sub_it = subDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
    {
      const SubGridEntityType& sub_entity = *sub_it;

      EntityPointerType host_entity_pointer = subGrid.template getHostEntity< 0 >(*sub_it);
      const EntityType& host_entity = *host_entity_pointer;

      SubGridLocalFunctionType sub_loc_value_1 = sub_func_1.localFunction(sub_entity);
      SubGridLocalFunctionType sub_loc_value_2 = sub_func_2.localFunction(sub_entity);
      LocalFunctionType host_loc_value_1 = host_func_1.localFunction(host_entity);
      LocalFunctionType host_loc_value_2 = host_func_2.localFunction(host_entity);

      const unsigned int numBaseFunctions = sub_loc_value_1.baseFunctionSet().size();
      for (unsigned int i = 0; i < numBaseFunctions; ++i)
      {
        host_loc_value_1[i] = sub_loc_value_1[i];
        host_loc_value_2[i] = sub_loc_value_2[i];
      }
    }
  } // subgrid_to_hostrid_function

  //! method to get the local mesh size H of a coarse grid entity 'T'
  // works only for our 2D examples!!!!
  RangeType get_coarse_grid_H(const EntityType& entity) {
    // entity_H means H (the diameter of the entity)
    RangeType entity_H = 0.0;

    const GridPartType& coarseGridPart = specifier_.coarseSpace().gridPart();

    // compute the size of the faces of the entities and selected the largest.
    IntersectionIteratorType endnit = coarseGridPart.iend(entity);

    for (IntersectionIteratorType nit = coarseGridPart.ibegin(entity); nit != endnit; ++nit)
    {
      FaceQuadratureType innerFaceQuadrature(coarseGridPart, *nit, 0, FaceQuadratureType::INSIDE);

      DomainType scaledOuterNormal
        = nit->integrationOuterNormal( innerFaceQuadrature.localPoint(0) );

      // get 'volume' of the visited face (this only works because we do not have curved faces):
      RangeType visitedFaceVolume(0.0);
      for (int k = 0; k < dimension; ++k)
        visitedFaceVolume += scaledOuterNormal[k] * scaledOuterNormal[k];
      visitedFaceVolume = sqrt(visitedFaceVolume);

      if (visitedFaceVolume > entity_H)
        entity_H = visitedFaceVolume;
    }

    return entity_H;
  } // get_coarse_grid_H

  // for a coarse grid entity T:
  // return:  H_T ||f||_{L^2(T)}
  RangeType indicator_f(const EntityType& entity) {
    // create quadrature for given geometry type
    CachingQuadrature< GridPartType, 0 > entityQuadrature(entity, 2 * spacePolOrd + 2);

    // get geoemetry of entity
    const EntityGeometryType& geometry = entity.geometry();

    RangeType H_T = get_coarse_grid_H(entity);

    RangeType y(0);
    RangeType local_indicator(0);

    const int quadratureNop = entityQuadrature.nop();
    for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
    {
      const double weight = entityQuadrature.weight(quadraturePoint)
                            * geometry.integrationElement( entityQuadrature.point(quadraturePoint) );

      f_.evaluate(geometry.global( entityQuadrature.point(quadraturePoint) ), y);
      y = y * y;

      local_indicator += weight * y;
    }

    local_indicator *= pow(H_T, 2.0);

    return sqrt(local_indicator);
  } // indicator_f

  // is a given point on a given face?
  bool point_on_face(const Intersection& face, const DomainType& point) {
    DomainType corner_0 = face.geometry().corner(0);
    DomainType corner_1 = face.geometry().corner(1);

    if (corner_0[0] == corner_1[0])
    {
      if (point[0] != corner_0[0])
      { return false; } else {
        RangeType lambda = (point[1] - corner_1[1]) / (corner_0[1] - corner_1[1]);
        if ( (lambda >= 0.0) && (lambda <= 1.0) )
        { return true; } else
        { return false; }
      }
    } else {
      RangeType lambda = (point[0] - corner_1[0]) / (corner_0[0] - corner_1[0]);

      if ( (lambda >= 0.0) && (lambda <= 1.0) )
      {
        RangeType convex_comb = (lambda * corner_0[1]) + ( (1.0 - lambda) * corner_1[1] );
        if (convex_comb == point[1])
        { return true; } else
        { return false; }
      } else
      { return false; }
    }
  } // point_on_face

  // is a given face part of another given coarse face
  bool is_subface(const Intersection& fine_face, const Intersection& coarse_face) {
    DomainType corner_0 = fine_face.geometry().corner(0);
    DomainType corner_1 = fine_face.geometry().corner(1);

    if ( point_on_face(coarse_face, corner_0) && point_on_face(coarse_face, corner_1) )
    { return true; } else
    { return false; }
  } // is_subface

  // jump in conservative flux
  void getFluxes(const EntityType& coarse_entity,
                 const DiscreteFunctionType& msfem_coarse_part,
                 RangeType& jump_conservative_flux,
                 RangeType& jump_coarse_flux) {
    // jump for each face
    RangeType jump[3];

    jump[0] = 0.0;
    jump[1] = 0.0;
    jump[2] = 0.0;

    // coarse grid jump for each face
    RangeType coarse_jump[3];

    coarse_jump[0] = 0.0;
    coarse_jump[1] = 0.0;
    coarse_jump[2] = 0.0;

    const DiscreteFunctionSpaceType& coarseDiscreteFunctionSpace = specifier_.coarseSpace();
    const LeafIndexSetType& coarseGridLeafIndexSet = coarseDiscreteFunctionSpace.gridPart().grid().leafIndexSet();

    int index_coarse_entity = coarseGridLeafIndexSet.index(coarse_entity);

    //! ---- get sub grid that for coarse entity

    // the sub grid U(T) that belongs to the coarse_grid_entity T
    SubGridType& sub_grid_U_T = subgrid_list_.getSubGrid(index_coarse_entity);
    SubGridPart subGridPart(sub_grid_U_T);

    SubGridDiscreteFunctionSpaceType localDiscreteFunctionSpace(subGridPart);

    SubGridDiscreteFunctionType conservative_flux_coarse_ent_e0("Conservative Flux on coarse entity for e_0",
                                                                localDiscreteFunctionSpace);
    conservative_flux_coarse_ent_e0.clear();

    SubGridDiscreteFunctionType conservative_flux_coarse_ent_e1("Conservative Flux on coarse entity for e_1",
                                                                localDiscreteFunctionSpace);
    conservative_flux_coarse_ent_e1.clear();

    // --------- load local solutions -------
    boost::format flux_location("%s/cf_problems/_conservativeFlux_e_%d_sg_%d");
    // the file/place, where we saved the solutions conservative flux problems problems
    const std::string cf_solution_location_0 = (flux_location
                                               % path_ % 0 % index_coarse_entity).str();
    // reader for data file:
    DiscreteFunctionReader(cf_solution_location_0).read(0, conservative_flux_coarse_ent_e0);
    // flux for e_1 ...
    const std::string cf_solution_location_1 = (flux_location
                                               % path_ % 1 % index_coarse_entity).str();
    // reader for data file:
    DiscreteFunctionReader(cf_solution_location_1).read(0, conservative_flux_coarse_ent_e1);

    DiscreteFunctionType cflux_coarse_ent_e0_host("Conservative Flux on coarse entity for e_0",
                                                  fineDiscreteFunctionSpace_);
    DiscreteFunctionType cflux_coarse_ent_e1_host("Conservative Flux on coarse entity for e_1",
                                                  fineDiscreteFunctionSpace_);

    subgrid_to_hostrid_function(conservative_flux_coarse_ent_e0,
                                conservative_flux_coarse_ent_e1,
                                cflux_coarse_ent_e0_host,
                                cflux_coarse_ent_e1_host);

    // flux for each neighbor entity
    //!TODO automatic memory
    DiscreteFunctionType* cflux_neighbor_ent_e0_host[3];
    DiscreteFunctionType* cflux_neighbor_ent_e1_host[3];

    std::vector< IntersectionIteratorType > coarse_face;
    RangeType coarse_face_volume[3];

    const GridPartType& coarseGridPart = specifier_.coarseSpace().gridPart();

    int local_face_index = 0;

    IntersectionIteratorType endnit = coarseGridPart.iend(coarse_entity);
    for (IntersectionIteratorType face_it = coarseGridPart.ibegin(coarse_entity); face_it != endnit; ++face_it)
    {
      coarse_face.push_back(face_it);

      coarse_face_volume[local_face_index] = face_it->geometry().volume();

      if ( face_it->neighbor() )
      {
        EntityPointerType outside_it = face_it->outside();

        int index_coarse_neighbor_entity = coarseGridLeafIndexSet.index(*outside_it);

        // --- get subgrids and load fluxes ---

        SubGridType& sub_grid_neighbor_U_T = subgrid_list_.getSubGrid(index_coarse_neighbor_entity);

        SubGridPart subGridPart_neighbor(sub_grid_neighbor_U_T);
        SubGridDiscreteFunctionSpaceType localDiscreteFunctionSpace_neighbor(subGridPart_neighbor);

        SubGridDiscreteFunctionType conservative_flux_coarse_ent_e0_neighbor(
          "Conservative Flux on neighbor coarse entity for e_0",
          localDiscreteFunctionSpace_neighbor);
        conservative_flux_coarse_ent_e0_neighbor.clear();

        SubGridDiscreteFunctionType conservative_flux_coarse_ent_e1_neighbor(
          "Conservative Flux on neighbor coarse entity for e_1",
          localDiscreteFunctionSpace_neighbor);
        conservative_flux_coarse_ent_e1_neighbor.clear();

        // --------- load local solutions -------
        // the file/place, where we saved the solutions conservative flux problems problems
        const std::string cf_solution_location_0_neighbor = (flux_location
                                          % path_ % 0 % index_coarse_neighbor_entity).str();
        // reader for data file:
        DiscreteFunctionReader(cf_solution_location_0_neighbor).read(0, conservative_flux_coarse_ent_e0_neighbor);

        // flux for e_1 ...
        const std::string cf_solution_location_1_neighbor = (flux_location
                                          % path_ % 1 % index_coarse_neighbor_entity).str();
        // reader for data file:
        DiscreteFunctionReader(cf_solution_location_1_neighbor).read(0, conservative_flux_coarse_ent_e1_neighbor);

        cflux_neighbor_ent_e0_host[local_face_index] = new DiscreteFunctionType(
          "Conservative Flux on neighbor coarse entity for e_0",
          fineDiscreteFunctionSpace_);

        cflux_neighbor_ent_e1_host[local_face_index] = new DiscreteFunctionType(
          "Conservative Flux on neighbor coarse entity for e_1",
          fineDiscreteFunctionSpace_);

        subgrid_to_hostrid_function(conservative_flux_coarse_ent_e0_neighbor,
                                    conservative_flux_coarse_ent_e1_neighbor,
                                    *cflux_neighbor_ent_e0_host[local_face_index],
                                    *cflux_neighbor_ent_e1_host[local_face_index]);
      }

      local_face_index += 1;
    }

    if (local_face_index != 3)
    {
      DUNE_THROW(Dune::InvalidStateException,"Error! Implementation only for triangular mesh in 2d!");
    }

    SubGridIteratorType sub_endit = localDiscreteFunctionSpace.end();
    for (SubGridIteratorType sub_it = localDiscreteFunctionSpace.begin(); sub_it != sub_endit; ++sub_it)
    {
      const SubGridEntityType& DUNE_UNUSED(sub_entity) = *sub_it;

      EntityPointerType host_entity_pointer = sub_grid_U_T.template getHostEntity< 0 >(*sub_it);
      const EntityType& host_entity = *host_entity_pointer;

      EntityPointerType father_of_sub_grid_entity = host_entity_pointer;
      for (int lev = 0; lev < specifier_.getLevelDifference(); ++lev)
        father_of_sub_grid_entity = father_of_sub_grid_entity->father();
      EntityPointerType coarse_father_test = father_of_sub_grid_entity;
      bool father_found = false;
      while (father_found == false)
      {
        if (coarseGridLeafIndexSet.contains(*coarse_father_test) == true)
        { father_of_sub_grid_entity = coarse_father_test; }

        if (coarse_father_test->hasFather() == false)
        { father_found = true; } else
        { coarse_father_test = coarse_father_test->father(); }
      }

      int coarse_sub_father_index = coarseGridLeafIndexSet.index(*father_of_sub_grid_entity);
      if (coarse_sub_father_index != index_coarse_entity)
      { continue; }

      LocalFunctionType loc_cf_coarse_ent_e0 = cflux_coarse_ent_e0_host.localFunction(host_entity);
      LocalFunctionType loc_cf_coarse_ent_e1 = cflux_coarse_ent_e1_host.localFunction(host_entity);

      LocalFunctionType loc_msfem_coarse_part = msfem_coarse_part.localFunction(host_entity);

      IntersectionIteratorType end_it_U_T = fineDiscreteFunctionSpace_.gridPart().iend(host_entity);
      for (IntersectionIteratorType face_it_U_T = fineDiscreteFunctionSpace_.gridPart().ibegin(host_entity);
           face_it_U_T != end_it_U_T; ++face_it_U_T)
      {
        int relevant_face_index = -1;

        if ( is_subface(*face_it_U_T, *coarse_face[0]) )
        { relevant_face_index = 0; }

        if ( is_subface(*face_it_U_T, *coarse_face[1]) )
        { relevant_face_index = 1; }

        if ( is_subface(*face_it_U_T, *coarse_face[2]) )
        { relevant_face_index = 2; }

        if ( (relevant_face_index == -1) || (face_it_U_T->neighbor() == false) )
        { continue; }

        EntityPointerType outside_sub_it = face_it_U_T->outside();

        LocalFunctionType loc_cf_coarse_neighbor_ent_e0
          = (*cflux_neighbor_ent_e0_host[relevant_face_index]).localFunction(host_entity);
        LocalFunctionType loc_cf_coarse_neighbor_ent_e1
          = (*cflux_neighbor_ent_e1_host[relevant_face_index]).localFunction(host_entity);

        LocalFunctionType loc_msfem_coarse_part_neighbor = msfem_coarse_part.localFunction(*outside_sub_it);

        // evaluate the gradient of the MsfEM coarse part in the center of the coarse entity
        EntityQuadratureType coarseEntQuadrature(host_entity, 0);
        JacobianRangeType gradient_msfem_coarse_ent(0.);
        loc_msfem_coarse_part.jacobian(coarseEntQuadrature[0], gradient_msfem_coarse_ent);

        // evaluate the gradient of the MsfEM coarse part in the center of the current neighbor of the coarse entity
        EntityQuadratureType coarseEntQuadratureNeighbor(*outside_sub_it, 0);
        JacobianRangeType gradient_msfem_coarse_neighbor_ent(0.);
        loc_msfem_coarse_part_neighbor.jacobian(coarseEntQuadratureNeighbor[0], gradient_msfem_coarse_neighbor_ent);

        FaceQuadratureType faceQuadrature(
          fineDiscreteFunctionSpace_.gridPart(), *face_it_U_T,
          2 * fineDiscreteFunctionSpace_.order() + 2, FaceQuadratureType::INSIDE);
        // inside macht hier keinen Unterschied, da wir formal stetige Funktionen haben und nicht die Gradienten
        // auswerten

        const FaceGeometryType& faceGeometry = face_it_U_T->geometry();

        const size_t numQuadraturePoints = faceQuadrature.nop();

        RangeType jump_integral(0.0);

        RangeType check_sum(0.0);
        for (size_t quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
        {
          typedef typename Intersection::LocalCoordinate LocalCoordinate;
          const LocalCoordinate local_point = faceGeometry.local( faceQuadrature.point(quadraturePoint) );

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

          JacobianRangeType coarse_jump_contribution;
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

    // std :: cout << "jump[0] = " << jump[0] << std :: endl;
    // std :: cout << "jump[1] = " << jump[1] << std :: endl;
    // std :: cout << "jump[2] = " << jump[2] << std :: endl;
    // std :: cout << std :: endl;

    jump_conservative_flux = ( sqrt(jump[0]) + sqrt(jump[1]) + sqrt(jump[2]) );
    jump_coarse_flux = sqrt(coarse_jump[0]) + sqrt(coarse_jump[1]) + sqrt(coarse_jump[2]);
  } // getFluxes

  // adaptive_refinement
  RangeType adaptive_refinement(GridType& /*coarse_grid*/,
                                const DiscreteFunctionType& msfem_solution,
                                const DiscreteFunctionType& msfem_coarse_part,
                                const DiscreteFunctionType& msfem_fine_part)
  {
    DSC_LOG_INFO << "Start computing conservative fluxes..." << std::endl;
    ConservativeFluxProblemSolver< SubGridDiscreteFunctionType, DiscreteFunctionType, DiffusionOperatorType,
                                   MacroMicroGridSpecifierImp >
    flux_problem_solver(fineDiscreteFunctionSpace_, diffusion_, specifier_, path_);

    flux_problem_solver.solve_all(subgrid_list_);

    DSC_LOG_INFO << "Conservative fluxes computed successfully." << std::endl;

    DSC_LOG_INFO << "Starting error estimation..." << std::endl;

    const int number_of_coarse_grid_entities = specifier_.getNumOfCoarseEntities();
    const int DUNE_UNUSED(number_of_fine_grid_entities) = msfem_solution.space().gridPart().grid().size(0 /*codim*/);

    // error contribution of H^2 ||f||_L2(T):
    std::vector< RangeType > loc_coarse_residual(number_of_coarse_grid_entities);

    std::vector< RangeType > loc_projection_error(number_of_coarse_grid_entities);

    std::vector< RangeType > loc_coarse_grid_jumps(number_of_coarse_grid_entities);

    std::vector< RangeType > loc_conservative_flux_jumps(number_of_coarse_grid_entities);

    // Summation ueber die einzelnen fine-grid indicators die zu coarse grid entity beitragen:

    std::vector< RangeType > loc_approximation_error(number_of_coarse_grid_entities);

    std::vector< RangeType > loc_fine_grid_jumps(number_of_coarse_grid_entities);

    specifier_.initialize_local_error_manager();

    for (int m = 0; m < number_of_coarse_grid_entities; ++m)
    {
      loc_coarse_residual[m] = 0.0;
      loc_projection_error[m] = 0.0;
      loc_coarse_grid_jumps[m] = 0.0;
      loc_conservative_flux_jumps[m] = 0.0;
      loc_approximation_error[m] = 0.0;
      loc_fine_grid_jumps[m] = 0.0;
    }

    RangeType total_coarse_residual(0.0);
    RangeType total_projection_error(0.0);
    RangeType total_coarse_grid_jumps(0.0);
    RangeType total_conservative_flux_jumps(0.0);
    RangeType total_approximation_error(0.0);
    RangeType total_fine_grid_jumps(0.0);

    RangeType total_estimated_error(0.0);

    const DiscreteFunctionSpaceType& coarseDiscreteFunctionSpace = specifier_.coarseSpace();
    const LeafIndexSetType& coarseGridLeafIndexSet = coarseDiscreteFunctionSpace.gridPart().grid().leafIndexSet();

    // Coarse Entity Iterator
    const IteratorType coarse_grid_end = coarseDiscreteFunctionSpace.end();
    for (IteratorType coarse_grid_it = coarseDiscreteFunctionSpace.begin();
         coarse_grid_it != coarse_grid_end;
         ++coarse_grid_it)
    {
      int global_index_entity = coarseGridLeafIndexSet.index(*coarse_grid_it);

      loc_coarse_residual[global_index_entity] = indicator_f(*coarse_grid_it);
      specifier_.set_loc_coarse_residual(global_index_entity, loc_coarse_residual[global_index_entity]);
      total_coarse_residual += pow(loc_coarse_residual[global_index_entity], 2.0);

      getFluxes(*coarse_grid_it,
                msfem_coarse_part,
                loc_conservative_flux_jumps[global_index_entity],
                loc_coarse_grid_jumps[global_index_entity]);

      specifier_.set_loc_coarse_grid_jumps(global_index_entity, loc_coarse_grid_jumps[global_index_entity]);

      specifier_.set_loc_conservative_flux_jumps(global_index_entity, loc_conservative_flux_jumps[global_index_entity]);

      total_conservative_flux_jumps += pow(loc_conservative_flux_jumps[global_index_entity], 2.0);
      total_coarse_grid_jumps += pow(loc_coarse_grid_jumps[global_index_entity], 2.0);

      // -------------------------- subgrids and local function ----------------------

      // the sub grid U(T) that belongs to the coarse_grid_entity T
      SubGridType& sub_grid_U_T = subgrid_list_.getSubGrid(global_index_entity);
      SubGridPart subGridPart(sub_grid_U_T);

      SubGridDiscreteFunctionSpaceType localDiscreteFunctionSpace(subGridPart);

      SubGridDiscreteFunctionType local_problem_solution_e0("Local problem Solution e_0", localDiscreteFunctionSpace);
      local_problem_solution_e0.clear();

      SubGridDiscreteFunctionType local_problem_solution_e1("Local problem Solution e_1", localDiscreteFunctionSpace);
      local_problem_solution_e1.clear();

      // --------- load local solutions -------
      // the file/place, where we saved the solutions of the cell problems
      const std::string local_solution_location = (boost::format("%s/local_problems/_localProblemSolutions_%d")
                                % path_ % global_index_entity).str();

      // reader for the cell problem data file:
      DiscreteFunctionReader discrete_function_reader(local_solution_location);
      discrete_function_reader.read(0, local_problem_solution_e0);
      discrete_function_reader.read(1, local_problem_solution_e1);

      // iterator for the local micro grid ('the subgrid corresponding with U(T)')
      const SubGridIteratorType local_grid_it_end = localDiscreteFunctionSpace.end();
      for (SubGridIteratorType local_grid_it = localDiscreteFunctionSpace.begin();
           local_grid_it != local_grid_it_end;
           ++local_grid_it)
      {
        const SubGridEntityType& local_grid_entity = *local_grid_it;

        // check if "local_grid_entity" (which is an entity of U(T)) is in T:
        // -------------------------------------------------------------------
        const EntityPointerType host_local_grid_it = localDiscreteFunctionSpace.grid().template getHostEntity< 0 >(
          local_grid_entity);

        EntityPointerType father_of_loc_grid_it = host_local_grid_it;

        for (int lev = 0; lev < specifier_.getLevelDifference(); ++lev)
          father_of_loc_grid_it = father_of_loc_grid_it->father();

        EntityPointerType coarse_father_test = father_of_loc_grid_it;
        bool father_found = false;
        while (father_found == false)
        {
          if (coarseGridLeafIndexSet.contains(*coarse_father_test) == true)
          { father_of_loc_grid_it = coarse_father_test; }

          if (coarse_father_test->hasFather() == false)
          { father_found = true; } else
          { coarse_father_test = coarse_father_test->father(); }
        }

        bool entities_identical = true;
        int number_of_nodes = (*coarse_grid_it).template count< 2 >();
        for (int k = 0; k < number_of_nodes; k += 1)
        {
          if ( !( coarse_grid_it->geometry().corner(k) == father_of_loc_grid_it->geometry().corner(k) ) )
          { entities_identical = false; }
        }

        if (entities_identical == false)
          continue;

        // -------------------------------------------------------------------
        EntityQuadratureType host_grid_quadrature(*host_local_grid_it, 0);
        LocalFunctionType localized_msfem_coarse_part = msfem_coarse_part.localFunction(*host_local_grid_it);
        LocalFunctionType localized_msfem_fine_part = msfem_fine_part.localFunction(*host_local_grid_it);

        JacobianRangeType grad_msfem_coarse_part;
        localized_msfem_coarse_part.jacobian(host_grid_quadrature[0], grad_msfem_coarse_part);

        JacobianRangeType grad_msfem_fine_part;
        localized_msfem_fine_part.jacobian(host_grid_quadrature[0], grad_msfem_fine_part);

        const SubGridEntityGeometryType& local_grid_geometry = local_grid_entity.geometry();
        assert(local_grid_entity.partitionType() == InteriorEntity);

        // quadrature formula on the sub grid entity

        // higher order quadrature, since A^{\epsilon} is highly variable
        LocalGridEntityQuadratureType local_grid_quadrature(local_grid_entity,
                                                            2 * localDiscreteFunctionSpace.order() + 2);
        const size_t numQuadraturePoints = local_grid_quadrature.nop();

        for (size_t localQuadraturePoint = 0; localQuadraturePoint < numQuadraturePoints; ++localQuadraturePoint)
        {
          // local (barycentric) coordinates (with respect to entity)
          const typename LocalGridEntityQuadratureType::CoordinateType& local_subgrid_point
            = local_grid_quadrature.point(localQuadraturePoint);

          DomainType global_point_in_U_T = local_grid_geometry.global(local_subgrid_point);

          const double weight_local_quadrature
            = local_grid_quadrature.weight(localQuadraturePoint) * local_grid_geometry.integrationElement(
            local_subgrid_point);

          SubGridLocalFunctionType localized_local_problem_solution_e0 = local_problem_solution_e0.localFunction(
            local_grid_entity);
          SubGridLocalFunctionType localized_local_problem_solution_e1 = local_problem_solution_e1.localFunction(
            local_grid_entity);

          // grad corrector for e_0 and e_1
          JacobianRangeType grad_loc_sol_e0, grad_loc_sol_e1;
          localized_local_problem_solution_e0.jacobian(local_grid_quadrature[localQuadraturePoint], grad_loc_sol_e0);
          localized_local_problem_solution_e1.jacobian(local_grid_quadrature[localQuadraturePoint], grad_loc_sol_e1);

          JacobianRangeType projection_error_gradient;
          for (int k = 0; k < dimension; ++k)
            projection_error_gradient[0][k] = grad_msfem_fine_part[0][k]
                                              - (grad_loc_sol_e0[0][k] * grad_msfem_coarse_part[0][0]
                                                 + grad_loc_sol_e1[0][k] * grad_msfem_coarse_part[0][1]);

          JacobianRangeType diffusive_flux_projection;
          diffusion_.diffusiveFlux(global_point_in_U_T, projection_error_gradient, diffusive_flux_projection);

          RangeType value(0.0);
          for (int k = 0; k < dimension; ++k)
          {
            value += pow(diffusive_flux_projection[0][k], 2.0);
          }

          loc_projection_error[global_index_entity] += value * weight_local_quadrature;
          total_projection_error += value * weight_local_quadrature;
        }
      }
    }

    // fine-grid iterator:

    // Coarse Entity Iterator
    const IteratorType fine_grid_end = fineDiscreteFunctionSpace_.end();
    for (IteratorType fine_grid_it = fineDiscreteFunctionSpace_.begin(); fine_grid_it != fine_grid_end; ++fine_grid_it)
    {
      EntityType& entity = *fine_grid_it;

      // identify coarse grid father entity
      EntityPointerType coarse_father(*fine_grid_it);
      for (int lev = 0; lev < specifier_.getLevelDifference(); ++lev)
        coarse_father = coarse_father->father();

      //! new version:
      EntityPointerType coarse_father_test = coarse_father;

      bool father_found = false;
      while (father_found == false)
      {
        if (coarseGridLeafIndexSet.contains(*coarse_father_test) == true)
        { coarse_father = coarse_father_test; }

        if (coarse_father_test->hasFather() == false)
        { father_found = true; } else
        { coarse_father_test = coarse_father_test->father(); }
      }

      int coarse_father_index = coarseGridLeafIndexSet.index(*coarse_father);

      const EntityGeometryType& entityGeometry = entity.geometry();

      EntityQuadratureType entityQuadrature(entity, 0);    // 0 = polynomial order
      const DomainType& x = entityGeometry.global( entityQuadrature.point(0) );

      LocalFunctionType local_msfem_sol = msfem_solution.localFunction(entity);
      JacobianRangeType gradient_msfem_sol(0.);
      local_msfem_sol.jacobian(entityQuadrature[0], gradient_msfem_sol);

      JacobianRangeType diffusive_flux_x;
      diffusion_.diffusiveFlux(x, gradient_msfem_sol, diffusive_flux_x);

      EntityQuadratureType highOrder_entityQuadrature(entity, 2 * spacePolOrd + 2);

      const int quadratureNop = highOrder_entityQuadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
      {
        const double weight = highOrder_entityQuadrature.weight(quadraturePoint)
                              * entityGeometry.integrationElement( highOrder_entityQuadrature.point(quadraturePoint) );

        DomainType point = entityGeometry.global( highOrder_entityQuadrature.point(quadraturePoint) );

        JacobianRangeType diffusive_flux_high_order;
        diffusion_.diffusiveFlux(point, gradient_msfem_sol, diffusive_flux_high_order);

        RangeType value = 0.0;
        for (int i = 0; i < dimension; ++i)
          value += pow(diffusive_flux_x[0][i] - diffusive_flux_high_order[0][i], 2.0);

        loc_approximation_error[coarse_father_index] += weight * value;
        total_approximation_error += weight * value;
      }

      const GridPartType& fineGridPart = fineDiscreteFunctionSpace_.gridPart();

      IntersectionIteratorType endnit = fineGridPart.iend(*fine_grid_it);
      for (IntersectionIteratorType nit = fineGridPart.ibegin(*fine_grid_it); nit != endnit; ++nit)
      {
        FaceQuadratureType innerFaceQuadrature(fineGridPart, *nit, 0, FaceQuadratureType::INSIDE);
        const FaceGeometryType& faceGeometry = nit->geometry();

        if (nit->neighbor() == false)
        { continue; }

        EntityPointerType outer_fine_grid_it = nit->outside();
        EntityType& outer_entity = *outer_fine_grid_it;

        EntityQuadratureType outer_entityQuadrature(outer_entity, 0);    // 0 = polynomial order
        const EntityGeometryType& outer_entityGeometry = outer_entity.geometry();
        const DomainType& outer_x = outer_entityGeometry.global( outer_entityQuadrature.point(0) );

        LocalFunctionType outer_local_msfem_sol = msfem_solution.localFunction(outer_entity);
        JacobianRangeType outer_gradient_msfem_sol(0.);
        outer_local_msfem_sol.jacobian(outer_entityQuadrature[0], outer_gradient_msfem_sol);

        JacobianRangeType diffusive_flux_outside;
        diffusion_.diffusiveFlux(outer_x, outer_gradient_msfem_sol, diffusive_flux_outside);

        DomainType unitOuterNormal
          = nit->unitOuterNormal( innerFaceQuadrature.localPoint(0) );

        const RangeType edge_length = faceGeometry.volume();

        RangeType int_value = 0.0;
        for (int i = 0; i < dimension; ++i)
          int_value[i] += (diffusive_flux_x[0][i] - diffusive_flux_outside[0][i]) * unitOuterNormal[i];
        int_value = pow(int_value, 2.0);

        loc_fine_grid_jumps[coarse_father_index] += edge_length * edge_length * int_value;
        total_fine_grid_jumps += edge_length * edge_length * int_value;
      }
    }

    for (int m = 0; m < number_of_coarse_grid_entities; ++m)
    {
      loc_approximation_error[m] = sqrt(loc_approximation_error[m]);
      loc_fine_grid_jumps[m] = sqrt(loc_fine_grid_jumps[m]);
      loc_projection_error[m] = sqrt(loc_projection_error[m]);

      specifier_.set_loc_approximation_error(m, loc_approximation_error[m]);
      specifier_.set_loc_fine_grid_jumps(m, loc_fine_grid_jumps[m]);
      specifier_.set_loc_projection_error(m, loc_projection_error[m]);
    }
    total_coarse_residual = sqrt(total_coarse_residual);
    total_projection_error = sqrt(total_projection_error);
    total_coarse_grid_jumps = sqrt(total_coarse_grid_jumps);
    total_conservative_flux_jumps = sqrt(total_conservative_flux_jumps);
    total_approximation_error = sqrt(total_approximation_error);
    total_fine_grid_jumps = sqrt(total_fine_grid_jumps);

    total_estimated_error += total_coarse_residual;
    total_estimated_error += total_projection_error;
    total_estimated_error += total_coarse_grid_jumps;
    total_estimated_error += total_conservative_flux_jumps;
    total_estimated_error += total_approximation_error;
    total_estimated_error += total_fine_grid_jumps;

    DSC_LOG_INFO  << std::endl
                  << "Estimated Errors:" << std::endl << std::endl
                  << "Total estimated error = " << total_estimated_error << "." << std::endl
                  << "where: " << std::endl
                  << "total_coarse_residual = " << total_coarse_residual << "." << std::endl
                  << "total_projection_error = " << total_projection_error << "." << std::endl
                  << "total_coarse_grid_jumps = " << total_coarse_grid_jumps << "." << std::endl
                  << "total_conservative_flux_jumps = " << total_conservative_flux_jumps << "." << std::endl
                  << "total_approximation_error = " << total_approximation_error << "." << std::endl
                  << "total_fine_grid_jumps = " << total_fine_grid_jumps << "." << std::endl;

    return total_estimated_error;
  } // adaptive_refinement
}; // end of class MsFEMErrorEstimator
} // end namespace

#endif // ifndef DUNE_MSFEM_ERRORESTIMATOR_HH
