// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MULTISCALE_IMPROVED_L2ERROR_HH
#define DUNE_MULTISCALE_IMPROVED_L2ERROR_HH

#include <config.h>
#include <dune/common/exceptions.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/multiscale/tools/misc/linear-lagrange-interpolation.hh>
#include <dune/multiscale/hmm/cell_problem_numbering.hh>
#include <dune/multiscale/tools/discretefunctionwriter.hh>

namespace Dune {

//! KLASSE IST NUR FUER DEN 2D FALL UND DIE VERWENDUNG VON EINER SIMPLIZIALEN TRIANGULIERUNG, SONST FUNKTIONIERT DAS
// ALLES NICHT!!!!!
#if 0
template< class DiscreteFunctionType, int n = 0 >
class ImprovedL2Error
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::RangeType                 RangeType;
  typedef typename DiscreteFunctionType::DomainType                DomainType;
  typedef typename DiscreteFunctionType::RangeFieldType            RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::IteratorType         IteratorType;

  typedef typename DiscreteFunctionSpaceType::GridType     GridType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  // typedef HierarchicIterator< GridType > HierarchicIteratorType;

// typename GridType :: HierarchicIteratorType    HierarchicIteratorType;

  typedef typename GridType::template Codim< 0 >::Entity        EntityType;
  typedef typename GridType::template Codim< 0 >::EntityPointer EntityPointerType;
  typedef typename GridType::template Codim< 0 >::Geometry      EnGeometryType;
  typedef typename GridType::Traits::GlobalIdSet                IdSet;

  typedef typename IdSet::IdType IdType;

  typedef typename EntityType::ctype coordType;

  // hierarchic iterator over all elements of codim 0
  typedef typename EntityType::HierarchicIterator               HierarchicIteratorType;
  typedef typename GridType::template Codim< 0 >::LevelIterator LevelIteratorType;

  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  enum { dim = GridType::dimension };
  enum { spacePolOrd = DiscreteFunctionSpaceType::polynomialOrder };
  enum { dimRange = RangeType::dimension };

public:
  //! for two discrete functions auf dem gleichen, aber nicht dem selben grid (z.B. identisch angelegt)
  template< int polOrd >
  RangeFieldType norm_L2(const DiscreteFunctionType& f1,
                         const DiscreteFunctionType& f2, double /*dummy*/ = 0) const {
    static_assert(dimRange==1,"L2Error::norm2: only implemented for dimRange = 1! \n");
    const DiscreteFunctionSpaceType& space = f1.space();
    const GridPartType& gridPart = space.gridPart();
    const auto& comm = gridPart.grid().comm();

    RangeFieldType ret = 0;
    // get function space
    const DiscreteFunctionSpaceType& dfsp = f1.space();
    const DiscreteFunctionSpaceType& dfsp_2 = f2.space();
    assert( dfsp.size() == dfsp_2.size() );

    RangeType lv1, lv2;
    // for product:
    // int quadOrd = (polOrd*dim)*(polOrd*dim);
    const int quadOrd = polOrd;

    // iterate over all elements defining the function
    const IteratorType eit = dfsp.end();
    IteratorType it_2 = dfsp_2.begin();
    for (IteratorType it = dfsp.begin(); it != eit; ++it, ++it_2)
    {
      const EntityType& en = *it;
      const EntityType& en_2 = *it_2;

      const Fem::CachingQuadrature< GridPartType, 0 > quad(en, quadOrd);
      // get local functions on current element
      const LocalFunctionType lf1 = f1.localFunction(en);
      const LocalFunctionType lf2 = f2.localFunction(en_2);

      // get geoemetry of entity
      const EnGeometryType& geo = en.geometry();

      const int quadNop = quad.nop();
      for (int qp = 0; qp < quadNop; ++qp)
      {
        const double det
          = geo.integrationElement( quad.point(qp) );
        // evaluate local functions
        lf1.evaluate(quad[qp], lv1);
        lf2.evaluate(quad[qp], lv2);
        // substract
        lv1 -= lv2;

        ret += det * quad.weight(qp) * (lv1 * lv1);
      } // end qp iteration
    } // end element iteration

    ret = comm.sum(ret);
    return sqrt(ret);
  } // end method

  // for the following method, the DiscreteFunction 'fine_disc_func' needs to be defined on a gridPart that is a
  // refinment of the gridPart of 'coarse_disc_func'
  template< int polOrd >
  RangeFieldType norm_uniform_grids(const DiscreteFunctionType& coarse_disc_func,
                                    const DiscreteFunctionType& fine_disc_func, double /*dummy*/ = 0) const {
    static_assert(dimRange==1,"L2Error::norm2: only implemented for dimRange = 1! \n");
    // get function spaces
    const DiscreteFunctionSpaceType& coarse_discreteFunctionSpace = coarse_disc_func.space();
    const DiscreteFunctionSpaceType& fine_discreteFunctionSpace = fine_disc_func.space();

    const GridPartType& fine_gridPart = fine_discreteFunctionSpace.gridPart();

    const IteratorType fine_element_reference = fine_discreteFunctionSpace.begin();
    const IteratorType coarse_element_reference = coarse_discreteFunctionSpace.begin();

    const int coarse_grid_level = coarse_element_reference->level();
    const int fine_grid_level = fine_element_reference->level();
    const int level_difference = fine_grid_level - coarse_grid_level;
    const auto& fine_comm = fine_gridPart.grid().comm();

    // to return the L2 Norm:
    RangeFieldType l2Norm = 0.0;
    // for product:
    // int quadOrd = (polOrd*2)*(polOrd*2);
    const int quadOrd = polOrd;

    // last entity of fine grid:
    const IteratorType fine_end = fine_discreteFunctionSpace.end();
    for (IteratorType fine_it = fine_discreteFunctionSpace.begin(); fine_it != fine_end; ++fine_it)
    {
      const EntityPointerType fine_father_entity = fine_it;
      for (int lev = 0; lev < level_difference; ++lev)
        fine_father_entity = fine_father_entity->father();

      const Fem::CachingQuadrature< GridPartType, 0 > fine_quad(*fine_it, quadOrd);

      // get local functions on current element
      const LocalFunctionType local_fine_disc_func = fine_disc_func.localFunction(*fine_it);
      const LocalFunctionType local_coarse_disc_func = coarse_disc_func.localFunction(*fine_father_entity);

      // create at quadrature with 3 quadrature points:
      const Fem::CachingQuadrature< GridPartType, 0 > coarse_quad(*fine_father_entity, 2);       // 3 points for linear pol in 2D
      // get geoemetry of coarse entity
      const EnGeometryType& coarse_geo = fine_father_entity->geometry();

      DomainType coarse_quad_point[3];
      RangeType local_value_coarse_func[3];
      for (int qp = 0; qp < 3; ++qp)
      {
        coarse_quad_point[qp] = coarse_geo.global( coarse_quad.point(qp) );

        // evaluate local coarse disc function
        local_coarse_disc_func.evaluate(coarse_quad[qp], local_value_coarse_func[qp]);
      }

      const LinearLagrangeFunction2D< DiscreteFunctionSpaceType >
      coarse_disc_func_entity(coarse_quad_point[0],
                              local_value_coarse_func[0],
                              coarse_quad_point[1],
                              local_value_coarse_func[1],
                              coarse_quad_point[2],
                              local_value_coarse_func[2]);

      // get geoemetry of entity
      const EnGeometryType& fine_geo = fine_it->geometry();

      const int fine_quadNop = fine_quad.nop();
      for (int qp = 0; qp < fine_quadNop; ++qp)
      {
        const double det
          = fine_geo.integrationElement( fine_quad.point(qp) );

        // find best quadrature point in the coarse grid quadrature:
        const DomainType fine_quad_point = fine_geo.global( fine_quad.point(qp) );

        RangeType coarse_value;
        coarse_disc_func_entity.evaluate(fine_quad_point, coarse_value);

        RangeType fine_value = 0.0;
        // evaluate fine local function
        local_fine_disc_func.evaluate(fine_quad[qp], fine_value);
        fine_value -= coarse_value;

        l2Norm += det * fine_quad.weight(qp) * (fine_value * fine_value);
      }       // end qp iteration
    }   // end fine grid element iteration

    l2Norm = fine_comm.sum(l2Norm);

    return sqrt(l2Norm);
  } // end method

  // expensive hack (no more required):
  template< int polOrd >
  RangeFieldType norm_adaptive_grids_2(const DiscreteFunctionType& coarse_disc_func,
                                       const DiscreteFunctionType& fine_disc_func, double /*dummy*/ = 0) const {
    // check if the discrete functions have valid dofs:
    if ( !coarse_disc_func.dofsValid() || !fine_disc_func.dofsValid() )
    {
      DSC_LOG_ERROR << "Solution of discrete function invalid." << std::endl;
      return 0.0;
    }

    bool error_in_compuation = false;

    // get function spaces
    const DiscreteFunctionSpaceType& coarse_discreteFunctionSpace = coarse_disc_func.space();
    const DiscreteFunctionSpaceType& fine_discreteFunctionSpace = fine_disc_func.space();
    const GridPartType& fine_gridPart = fine_discreteFunctionSpace.gridPart();
    const auto& fine_comm = fine_gridPart.grid().comm();

    // to return the L2 Norm:
    RangeFieldType l2Norm = 0.0;

    if (dimRange > 1)
    {
      DUNE_THROW(Dune::NotImplemented,"L2Error::norm2: only implemented for dimRange = 1! \n");
    }

    // das Zentrum des Referenzelements (bei einer Triangulierung mit Dreiecken): (1/3 , 1/3)
    DomainType center_of_reference_element;
    center_of_reference_element[0] = (1.0 / 3.0);
    center_of_reference_element[1] = (1.0 / 3.0);

    // die Eckpunkte des Referenzelements ( (0,0), (1,0) und (0,1) )

    // (0,0)
    DomainType reference_corner_0;
    reference_corner_0[0] = 0.0;
    reference_corner_0[1] = 0.0;

    // (0,1)
    DomainType reference_corner_1;
    reference_corner_1[0] = 0.0;
    reference_corner_1[1] = 1.0;
    // map to global

    // (1,0)
    DomainType reference_corner_2;
    reference_corner_2[0] = 1.0;
    reference_corner_2[1] = 0.0;

    // for product:
    // int quadOrd = (polOrd*2)*(polOrd*2);
    int quadOrd = polOrd;

    // last entity of fine grid:
    const IteratorType fine_end = fine_discreteFunctionSpace.end();
    for (IteratorType fine_it = fine_discreteFunctionSpace.begin(); fine_it != fine_end; ++fine_it)
    {
      // get geoemetry of fine grid entity:
      const EnGeometryType& fine_geo = fine_it->geometry();

      // das Zentrum des Fine-Grid Elements:
      const DomainType center_of_fine_it = fine_geo.global(center_of_reference_element);

      // wir klappern jetzt alle Makroelemente ab, um das zum fine-grid-element gehörende, relevante coarse-grid
      // element zu finden
      // dort werden die passenden Werte fuer die Lagrange-Interpolation bestimmt.

      // wurde das relevante coarse-grid element schon gefunden?
      bool relevant_coarse_entity_found = false;

      // for the Lagrange Interpolation on the relevant coarse entity:
      DomainType coarse_quad_point[3];
      RangeType local_value_coarse_func[3];

      const IteratorType coarse_entity_end = coarse_discreteFunctionSpace.end();
      for (IteratorType coarse_it = coarse_discreteFunctionSpace.begin(); coarse_it != coarse_entity_end; ++coarse_it)
      {
        // get geoemetry of coarse entity
        const EnGeometryType& coarse_geo = coarse_it->geometry();

        const LocalFunctionType local_coarse_disc_func = coarse_disc_func.localFunction(*coarse_it);

        // beschreibe das Dreieck mit Hilfe seiner Eckpunkt-Koordinaten (Konvexkombination) und schaue ob das Zentrum
        // des Fine-Grid elements in dieser Konvexkombination liegt.
        // map the reference corners to the global corners:
        //! TODO static check for simplicial grid
        const DomainType global_corner_0 = coarse_geo.global(reference_corner_0);
        const DomainType global_corner_1 = coarse_geo.global(reference_corner_1);
        const DomainType global_corner_2 = coarse_geo.global(reference_corner_2);

        // sei c (center) der Punkt von dem wir testen wollen, ob er im Coarse Grid Dreieck liegt, dann muss gelten
        // (Konvexkombination), dass
        // lambda_0 und lambda_1 existieren, mit (0<=lambda_0<=1) && (0<=lambda_1<=1) && (0<=lambda_0+lambda_1<=1) mit
        // c = e_2 + lambda_0 (e_0 - e_2) + lambda_1 (e_1 - e_2)
        // wobei e_0, e_1 und e_2 die Ecken des Dreiecks sind.

        RangeType lambda_0, lambda_1;

        // stellen wir das Gleichungssystem nach (lambda_0,lambda_1) auf, so ergibt sich:

        if ( (global_corner_0[0] - global_corner_2[0]) == 0 )
        {
          lambda_1 = (center_of_fine_it[0] - global_corner_2[0]) / (global_corner_1[0] - global_corner_2[0]);

          lambda_0
            = ( (center_of_fine_it[1] - global_corner_2[1]) + ( lambda_1 * (global_corner_2[1] - global_corner_1[1]) ) )
              / (global_corner_0[1] - global_corner_2[1]);
        } else {
          if ( (global_corner_1[0] - global_corner_2[0]) == 0 )
          {
            lambda_0 = (center_of_fine_it[0] - global_corner_2[0]) / (global_corner_0[0] - global_corner_2[0]);

            lambda_1
              = ( (center_of_fine_it[1] - global_corner_2[1]) - ( lambda_0 * (global_corner_0[1] - global_corner_2[1]) ) )
                / (global_corner_1[1] - global_corner_2[1]);
          } else {
            if ( (global_corner_1[1] - global_corner_2[1]) == 0 )
            {
              lambda_0 = (center_of_fine_it[1] - global_corner_2[1]) / (global_corner_0[1] - global_corner_2[1]);

              lambda_1
                = ( (center_of_fine_it[0] - global_corner_2[0]) - ( lambda_0 * (global_corner_0[0] - global_corner_2[0]) ) )
                  / (global_corner_1[0] - global_corner_2[0]);
            } else {
              lambda_1 = ( (center_of_fine_it[1] - global_corner_2[1]) / (global_corner_1[1] - global_corner_2[1]) )
                         - ( ( (global_corner_0[1] - global_corner_2[1]) / (global_corner_0[0] - global_corner_2[0]) )
                             * ( (center_of_fine_it[0] - global_corner_2[0]) / (global_corner_1[1] - global_corner_2[1]) ) );

              lambda_1 = lambda_1
                         / ( 1.0
                             - ( ( (global_corner_0[1] - global_corner_2[1]) / (global_corner_0[0] - global_corner_2[0]) )
                                 * ( (global_corner_1[0]
                                      - global_corner_2[0]) / (global_corner_1[1] - global_corner_2[1]) ) ) );

              lambda_0
                = ( (center_of_fine_it[0] - global_corner_2[0]) - ( lambda_1 * (global_corner_1[0] - global_corner_2[0]) ) )
                  / (global_corner_0[0] - global_corner_2[0]);
            }
          }
        }

        if ( (0.0 <= lambda_0) && (lambda_0 <= 1.0)
             && (0.0 <= lambda_1) && (lambda_1 <= 1.0)
             && ( (lambda_0 + lambda_1) <= 1.0 ) )
        {
          relevant_coarse_entity_found = true;

          // create at quadrature with 3 quadrature points:
          Fem::CachingQuadrature< GridPartType, 0 > coarse_quad(*coarse_it, 2);        // 3 points for linear pol in 2D
          for (int qp = 0; qp < 3; ++qp)
          {
            coarse_quad_point[qp] = coarse_geo.global( coarse_quad.point(qp) );

            // evaluate local coarse disc function
            local_coarse_disc_func.evaluate(coarse_quad[qp], local_value_coarse_func[qp]);
          }

          break;
        }
      }

      if (!relevant_coarse_entity_found)
      {
        DSC_LOG_ERROR << "In class >>ImprovedL2Error<<, in method >>norm_adaptive_grids<< :" << std::endl
                      << "No corresponding coarse grid entity found for fine grid entity => Error in computation of L2 error."
                      << std::endl;
        // DSC_LOG_ERROR << "Problem with fine-grid center: center_of_fine_it(" << center_of_fine_it[0] << "," <<
        // center_of_fine_it[1] << ")" << std :: endl;
        // if ( (center_of_fine_it[0]>= 0) || (center_of_fine_it[1]>= 0) )
        //! TODO throw??!?
        error_in_compuation = true;
      }

      //! determine the Lagrange Interpolation of the coarse discrete function on the relevant coarse-grid element
      // (determined in the previous loop):

      LinearLagrangeFunction2D< DiscreteFunctionSpaceType >
      coarse_disc_func_entity(coarse_quad_point[0],
                              local_value_coarse_func[0],
                              coarse_quad_point[1],
                              local_value_coarse_func[1],
                              coarse_quad_point[2],
                              local_value_coarse_func[2]);

      //! the fine grid quadrature:

      const Fem::CachingQuadrature< GridPartType, 0 > fine_quad(*fine_it, quadOrd);

      // get local functions on current element
      const LocalFunctionType local_fine_disc_func = fine_disc_func.localFunction(*fine_it);

      const int fine_quadNop = fine_quad.nop();
      for (int qp = 0; qp < fine_quadNop; ++qp)
      {
        const double det
          = fine_geo.integrationElement( fine_quad.point(qp) );

        // find best quadrature point in the coarse grid quadrature:
        const DomainType fine_quad_point = fine_geo.global( fine_quad.point(qp) );

        RangeType coarse_value;
        coarse_disc_func_entity.evaluate(fine_quad_point, coarse_value);

        RangeType fine_value = 0.0;
        // evaluate fine local function
        local_fine_disc_func.evaluate(fine_quad[qp], fine_value);
        fine_value -= coarse_value;

        l2Norm += det * fine_quad.weight(qp) * (fine_value * fine_value);
      }       // end qp iteration
    }   // end fine grid element iteration
    l2Norm = fine_comm.sum(l2Norm);

    if (error_in_compuation)
    {
      return 0.0;
    }
    return sqrt(l2Norm);
  } // end method


  // expensive hack
  // does not yet work:
  // evaluate  disc_func + Q^eps( disc_func ) in x
  template< typename LocalProblemNumManagerType, int polOrd >
  RangeType evaluate_with_corrector(const DiscreteFunctionType& coarse_disc_func,
                                    const DomainType& x,
                                    const DiscreteFunctionSpaceType& local_discreteFunctionSpace,
                                    const LocalProblemNumManagerType& lp_num_manager,
                                    double /*dummy*/ = 0) const {
    // f(x) + Q(f)(x)
    RangeType value = 0.0;

    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader(lp_num_manager.get_location());

    bool error_in_compuation = false;
    // get function spaces
    const DiscreteFunctionSpaceType& coarse_discreteFunctionSpace = coarse_disc_func.space();

    // die Eckpunkte des Referenzelements ( (0,0), (1,0) und (0,1) )
    // (0,0)
    DomainType reference_corner_0;
    reference_corner_0[0] = 0.0;
    reference_corner_0[1] = 0.0;

    // (1,0)
    DomainType reference_corner_1;
    reference_corner_1[0] = 1.0;
    reference_corner_1[1] = 0.0;
    // map to global

    // (0,1)
    DomainType reference_corner_2;
    reference_corner_2[0] = 0.0;
    reference_corner_2[1] = 1.0;

    // last entity of coarse grid:
    const IteratorType coarse_end = coarse_discreteFunctionSpace.end();
    for (IteratorType coarse_it = coarse_discreteFunctionSpace.begin(); coarse_it != coarse_end; ++coarse_it)
    {
      // T = coarse_it

      // get geoemetry of coarse grid entity:
      const EnGeometryType& geometry_of_T = coarse_it->geometry();

      // corner of the global element T:
      // ( map the reference corners to the global corners: )
      const DomainType corner_0_of_T = geometry_of_T.global(reference_corner_0);
      const DomainType corner_1_of_T = geometry_of_T.global(reference_corner_1);
      const DomainType corner_2_of_T = geometry_of_T.global(reference_corner_2);

      // Wir bilden die Referenzabbildung F : T_0 -> T, gegeben durch
      // F( y ) = A y + b
      // und überprüfen ob F^(-1)(x) in T_0 liegt

      // value of the matrix A (in F(x) = Ax + a_0) (belonging to the element T)
      double val_A_T[2][2];
      val_A_T[0][0] = (corner_1_of_T[0] - corner_0_of_T[0]);
      val_A_T[0][1] = (corner_2_of_T[0] - corner_0_of_T[0]);
      val_A_T[1][0] = (corner_1_of_T[1] - corner_0_of_T[1]);
      val_A_T[1][1] = (corner_2_of_T[1] - corner_0_of_T[1]);

      // define 'c := (a_1(1) - a_0(1))·(a_2(2) - a_0(2)) - (a_1(2) - a_0(2))·(a_2(1) - a_0(1))
      double c = 1.0 / ( (val_A_T[0][0] * val_A_T[1][1]) - (val_A_T[0][1] * val_A_T[1][0]) );
      // then the inverse A^{-1} is given by:
      // A^{-1}_11 = (1/c) (a_2(2) - a_0(2))     A^{-1}_12 = (1/c) (a_0(1) - a_2(1))
      // A^{-1}_21 = (1/c) (a_0(2) - a_1(2))     A^{-1}_22 = (1/c) (a_1(1) - a_0(1))

      // (A^{-1})^T:
      double val_A_T_inverse[2][2];
      val_A_T_inverse[0][0] = c * val_A_T[1][1];
      val_A_T_inverse[1][0] = c * (-1.0) * val_A_T[1][0];
      val_A_T_inverse[0][1] = c * (-1.0) * val_A_T[0][1];
      val_A_T_inverse[1][1] = c * val_A_T[0][0];

      // x = F( y ) = Ay + a_0
      // =>
      // y = A^{-1} ( x - a_0 )
      DomainType F_inverse_of_x(0.0);
      for (int k = 0; k < 2; ++k)
        for (int l = 0; l < 2; ++l)
          F_inverse_of_x[k] += val_A_T_inverse[k][l] * (x[l] - corner_0_of_T[l]);

      // y \in T_0 ?
      if ( ( (F_inverse_of_x[0] >= 0.0) && (F_inverse_of_x[0] <= 1.0) )
           && ( (F_inverse_of_x[1] >= 0.0) && ( F_inverse_of_x[1] <= (1.0 - F_inverse_of_x[0]) ) ) )
      {
        const LocalFunctionType local_coarse_disc_func = coarse_disc_func.localFunction(*coarse_it);

        // beschreibe das Dreieck mit Hilfe seiner Eckpunkt-Koordinaten (Konvexkombination) und schaue ob das Zentrum
        // des Fine-Grid elements in dieser Konvexkombination liegt.

        //! determine the Lagrange Interpolation of the coarse discrete function on the relevant coarse-grid element
        // (determined in the previous loop):

        // create at quadrature with 3 quadrature points:
        const Fem::CachingQuadrature< GridPartType, 0 > coarse_quad(*coarse_it, 2);       // 3 points for linear pol in 2D
        // for the Lagrange Interpolation on the relevant coarse entity:
        DomainType coarse_quad_point[3];
        RangeType local_value_coarse_func[3];

        for (int qp = 0; qp < 3; ++qp)
        {
          coarse_quad_point[qp] = geometry_of_T.global( coarse_quad.point(qp) );

          // evaluate local coarse disc function
          local_coarse_disc_func.evaluate(coarse_quad[qp], local_value_coarse_func[qp]);
        }

        LinearLagrangeFunction2D< DiscreteFunctionSpaceType >
        coarse_disc_func_entity(coarse_quad_point[0],
                                local_value_coarse_func[0],
                                coarse_quad_point[1],
                                local_value_coarse_func[1],
                                coarse_quad_point[2],
                                local_value_coarse_func[2]);

        coarse_disc_func_entity.evaluate(x, value);

        // corrector of the coarse function on the macro element T:
        DiscreteFunctionType corrector_of_coarse_disc_func("Corrector of the coarse function",
                                                           local_discreteFunctionSpace);
        corrector_of_coarse_disc_func.clear();

        DiscreteFunctionType corrector_of_base_func("Corrector of macro base function", local_discreteFunctionSpace);
        corrector_of_base_func.clear();

        typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
        const BasisFunctionSetType& baseSet = coarse_discreteFunctionSpace.basisFunctionSet(*coarse_it);
        const unsigned int numMacroBaseFunctions = baseSet.numBaseFunctions();

        int cell_problem_id[numMacroBaseFunctions];
        for (unsigned int i = 0; i < numMacroBaseFunctions; ++i)
        {
          cell_problem_id[i] = lp_num_manager.get_number_of_local_problem(coarse_it, i);
          discrete_function_reader.read(cell_problem_id[i], corrector_of_base_func);
          corrector_of_base_func *= local_coarse_disc_func[i];
          corrector_of_coarse_disc_func += corrector_of_base_func;
          corrector_of_base_func.clear();
        }

        const IteratorType local_grid_end = local_discreteFunctionSpace.end();
        for (IteratorType local_grid_it = local_discreteFunctionSpace.begin();
             local_grid_it != local_grid_end;
             ++local_grid_it)
        {
          // Let 'T' denote the coarse grid element in \Omega
          // with mapping F_T : T_0 -> T
          // and
          // let 'S' denote the local grid element in T_0 ( = F_T^{-1}(T_0) )
          // with mapping G_S : T_0 -> S
          // Then the global position (in \Omega) of a point y \in T_0 is given by
          // z = F_T ( y )

          // now, for 'x' we need to verify, if it is in F_T(S)
          // so we want to check:   F_T^{-1}(x) \in S
          // or alternatively:      G_S^{-1}( F_T^{-1}(x) ) \in G_S^{-1}(S) = T_0

          // we summarize:
          // check if:      G_S^{-1}( F_T^{-1}(x) ) \in T_0

          // S = local_grid_it

          // get geoemetry of the local grid entity:
          const EnGeometryType& geometry_of_S = local_grid_it->geometry();

          // corner of the local element S:
          // ( map the reference corners to the global corners: )
          const DomainType corner_0_of_S = geometry_of_S.global(reference_corner_0);
          const DomainType corner_1_of_S = geometry_of_S.global(reference_corner_1);
          const DomainType corner_2_of_S = geometry_of_S.global(reference_corner_2);

          // Wir bilden die Referenzabbildung G_S : T_0 -> S, gegeben durch
          // G_S( y ) = A_S y + b_S
          // und überprüfen ob G_S^(-1)( F_T^{-1}(x) ) in T_0 liegt

          // value of the matrix A_S (in G(x) = A_S x + a_0) (belonging to the element S)
          double val_A_S[2][2];
          val_A_S[0][0] = corner_1_of_S[0] - corner_0_of_S[0];
          val_A_S[0][1] = corner_2_of_S[0] - corner_0_of_S[0];
          val_A_S[1][0] = corner_1_of_S[1] - corner_0_of_S[1];
          val_A_S[1][1] = corner_2_of_S[1] - corner_0_of_S[1];

          // define 'c := (a_1(1) - a_0(1))·(a_2(2) - a_0(2)) - (a_1(2) - a_0(2))·(a_2(1) - a_0(1))
          const double c_S = 1.0 / ( (val_A_S[0][0] * val_A_S[1][1]) - (val_A_S[0][1] * val_A_S[1][0]) );
          // then the inverse A^{-1} is given by:
          // A^{-1}_11 = (1/c) (a_2(2) - a_0(2))     A^{-1}_12 = (1/c) (a_0(1) - a_2(1))
          // A^{-1}_21 = (1/c) (a_0(2) - a_1(2))     A^{-1}_22 = (1/c) (a_1(1) - a_0(1))

          // (A_S^{-1})^T:
          double val_A_S_inverse[2][2];
          val_A_S_inverse[0][0] = c_S * val_A_S[1][1];
          val_A_S_inverse[1][0] = c_S * (-1.0) * val_A_S[1][0];
          val_A_S_inverse[0][1] = c_S * (-1.0) * val_A_S[0][1];
          val_A_S_inverse[1][1] = c_S * val_A_S[0][0];

          // x = F_T( G_S(y) )
          // =>
          // F_T^{-1}(x) = G_S(y) = A_S y + (a_S)_0
          // =>
          // y = A_S^{-1} ( F_T^{-1}(x) - (a_S)_0 )

          // of this 'y', we need to check, if it belongs to T_0

          DomainType G_inverse_of_F_inverse_of_x(0.0);
          for (int k = 0; k < 2; ++k)
            for (int l = 0; l < 2; ++l)
              G_inverse_of_F_inverse_of_x[k] += val_A_S_inverse[k][l] * (F_inverse_of_x[l] - corner_0_of_S[l]);

          // y \in T_0 ?
          if ( ( (G_inverse_of_F_inverse_of_x[0] >= 0.0) && (G_inverse_of_F_inverse_of_x[0] <= 1.0) )
               && ( (G_inverse_of_F_inverse_of_x[1] >= 0.0)
                    && ( G_inverse_of_F_inverse_of_x[1] <= (1.0 - G_inverse_of_F_inverse_of_x[0]) ) ) )
          {
            LocalFunctionType local_corrector_disc_func = corrector_of_coarse_disc_func.localFunction(*local_grid_it);

            // create at quadrature with 3 quadrature points:
            Fem::CachingQuadrature< GridPartType, 0 > local_quad(*local_grid_it, 2);          // 3 points for linear pol in
                                                                                         // 2D
            // for the Lagrange Interpolation on the relevant coarse entity:
            DomainType fine_quad_point[3];
            RangeType local_value_corrector_func[3];

            for (int qp = 0; qp < 3; ++qp)
            {
              fine_quad_point[qp] = geometry_of_S.global( local_quad.point(qp) );

              // evaluate local coarse disc function
              local_corrector_disc_func.evaluate(local_quad[qp], local_value_corrector_func[qp]);
            }

            LinearLagrangeFunction2D< DiscreteFunctionSpaceType >
            corrector_disc_func_entity(fine_quad_point[0],
                                       local_value_corrector_func[0],
                                       fine_quad_point[1],
                                       local_value_corrector_func[1],
                                       fine_quad_point[2],
                                       local_value_corrector_func[2]);

            RangeType corrector_value = 0.0;
            corrector_disc_func_entity.evaluate(x, corrector_value);
            value += corrector_value;

            break;
          }
        }

        break;
      }
    }  // end coarse grid element iteration

    if (error_in_compuation)
    {
      return 0.0;
    }

    return value;
  } // end method

  // expensive hack
  // does not yet work:
  // evaluate gradient( disc_func + Q^eps( disc_func ) ) in x
  template< typename LocalProblemNumManagerType, int polOrd >
  RangeType jacobian_with_corrector(const DiscreteFunctionType& coarse_disc_func,
                                    const DomainType& x,
                                    const DiscreteFunctionSpaceType& local_discreteFunctionSpace,
                                    const LocalProblemNumManagerType& lp_num_manager,
                                    double /*dummy*/ = 0) const {
    DUNE_THROW(Dune::NotImplemented, "jacobian_with_corrector hack not working yet");
    // f(x) + Q(f)(x)
    RangeType value = 0.0;

    // reader for the cell problem data file:
    DiscreteFunctionReader discrete_function_reader(lp_num_manager.get_location());

    bool error_in_compuation = false;

    // get function spaces
    const DiscreteFunctionSpaceType& coarse_discreteFunctionSpace = coarse_disc_func.space();

    // die Eckpunkte des Referenzelements ( (0,0), (1,0) und (0,1) )
    // (0,0)
    DomainType reference_corner_0;
    reference_corner_0[0] = 0.0;
    reference_corner_0[1] = 0.0;

    // (1,0)
    DomainType reference_corner_1;
    reference_corner_1[0] = 1.0;
    reference_corner_1[1] = 0.0;
    // map to global

    // (0,1)
    DomainType reference_corner_2;
    reference_corner_2[0] = 0.0;
    reference_corner_2[1] = 1.0;

    // last entity of coarse grid:
    IteratorType coarse_end = coarse_discreteFunctionSpace.end();
    for (IteratorType coarse_it = coarse_discreteFunctionSpace.begin(); coarse_it != coarse_end; ++coarse_it)
    {
      // T = coarse_it

      // get geoemetry of coarse grid entity:
      const EnGeometryType& geometry_of_T = coarse_it->geometry();

      // corner of the global element T:
      // ( map the reference corners to the global corners: )
      const DomainType corner_0_of_T = geometry_of_T.global(reference_corner_0);
      const DomainType corner_1_of_T = geometry_of_T.global(reference_corner_1);
      const DomainType corner_2_of_T = geometry_of_T.global(reference_corner_2);

      // Wir bilden die Referenzabbildung F : T_0 -> T, gegeben durch
      // F( y ) = A y + b
      // und überprüfen ob F^(-1)(x) in T_0 liegt

      // value of the matrix A (in F(x) = Ax + a_0) (belonging to the element T)
      double val_A_T[2][2];
      val_A_T[0][0] = (corner_1_of_T[0] - corner_0_of_T[0]);
      val_A_T[0][1] = (corner_2_of_T[0] - corner_0_of_T[0]);
      val_A_T[1][0] = (corner_1_of_T[1] - corner_0_of_T[1]);
      val_A_T[1][1] = (corner_2_of_T[1] - corner_0_of_T[1]);

      // define 'c := (a_1(1) - a_0(1))·(a_2(2) - a_0(2)) - (a_1(2) - a_0(2))·(a_2(1) - a_0(1))
      const double c = 1.0 / ( (val_A_T[0][0] * val_A_T[1][1]) - (val_A_T[0][1] * val_A_T[1][0]) );
      // then the inverse A^{-1} is given by:
      // A^{-1}_11 = (1/c) (a_2(2) - a_0(2))     A^{-1}_12 = (1/c) (a_0(1) - a_2(1))
      // A^{-1}_21 = (1/c) (a_0(2) - a_1(2))     A^{-1}_22 = (1/c) (a_1(1) - a_0(1))

      // (A^{-1})^T:
      double val_A_T_inverse[2][2];
      val_A_T_inverse[0][0] = c * val_A_T[1][1];
      val_A_T_inverse[1][0] = c * (-1.0) * val_A_T[1][0];
      val_A_T_inverse[0][1] = c * (-1.0) * val_A_T[0][1];
      val_A_T_inverse[1][1] = c * val_A_T[0][0];

      // x = F( y ) = Ay + a_0
      // =>
      // y = A^{-1} ( x - a_0 )
      DomainType F_inverse_of_x(0.0);
      for (int k = 0; k < 2; ++k)
        for (int l = 0; l < 2; ++l)
          F_inverse_of_x[k] += val_A_T_inverse[k][l] * (x[l] - corner_0_of_T[l]);

      // y \in T_0 ?
      if ( ( (F_inverse_of_x[0] >= 0.0) && (F_inverse_of_x[0] <= 1.0) )
           && ( (F_inverse_of_x[1] >= 0.0) && ( F_inverse_of_x[1] <= (1.0 - F_inverse_of_x[0]) ) ) )
      {

        const LocalFunctionType local_coarse_disc_func = coarse_disc_func.localFunction(*coarse_it);

        // beschreibe das Dreieck mit Hilfe seiner Eckpunkt-Koordinaten (Konvexkombination) und schaue ob das Zentrum
        // des Fine-Grid elements in dieser Konvexkombination liegt.

        //! determine the Lagrange Interpolation of the coarse discrete function on the relevant coarse-grid element
        // (determined in the previous loop):

        // create at quadrature with 3 quadrature points:
        const Fem::CachingQuadrature< GridPartType, 0 > coarse_quad(*coarse_it, 2);       // 3 points for linear pol in 2D
        // for the Lagrange Interpolation on the relevant coarse entity:
        DomainType coarse_quad_point[3];
        RangeType local_value_coarse_func[3];

        for (int qp = 0; qp < 3; ++qp)
        {
          coarse_quad_point[qp] = geometry_of_T.global( coarse_quad.point(qp) );

          // evaluate local coarse disc function
          local_coarse_disc_func.evaluate(coarse_quad[qp], local_value_coarse_func[qp]);
        }

        LinearLagrangeFunction2D< DiscreteFunctionSpaceType >
        coarse_disc_func_entity(coarse_quad_point[0],
                                local_value_coarse_func[0],
                                coarse_quad_point[1],
                                local_value_coarse_func[1],
                                coarse_quad_point[2],
                                local_value_coarse_func[2]);

        coarse_disc_func_entity.evaluate(x, value);

        // corrector of the coarse function on the macro element T:
        DiscreteFunctionType corrector_of_coarse_disc_func("Corrector of the coarse function",
                                                           local_discreteFunctionSpace);
        corrector_of_coarse_disc_func.clear();

        DiscreteFunctionType corrector_of_base_func("Corrector of macro base function", local_discreteFunctionSpace);
        corrector_of_base_func.clear();

        typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
        const BasisFunctionSetType& baseSet = coarse_discreteFunctionSpace.basisFunctionSet(*coarse_it);
        const unsigned int numMacroBaseFunctions = baseSet.numBaseFunctions();

        int cell_problem_id[numMacroBaseFunctions];
        for (unsigned int i = 0; i < numMacroBaseFunctions; ++i)
        {
          cell_problem_id[i] = lp_num_manager.get_number_of_local_problem(coarse_it, i);

          discrete_function_reader.read(cell_problem_id[i], corrector_of_base_func);

          corrector_of_base_func *= local_coarse_disc_func[i];

          corrector_of_coarse_disc_func += corrector_of_base_func;

          corrector_of_base_func.clear();
        }

        IteratorType local_grid_end = local_discreteFunctionSpace.end();
        for (IteratorType local_grid_it = local_discreteFunctionSpace.begin();
             local_grid_it != local_grid_end;
             ++local_grid_it)
        {
          // Let 'T' denote the coarse grid element in \Omega
          // with mapping F_T : T_0 -> T
          // and
          // let 'S' denote the local grid element in T_0 ( = F_T^{-1}(T_0) )
          // with mapping G_S : T_0 -> S
          // Then the global position (in \Omega) of a point y \in T_0 is given by
          // z = F_T ( y )

          // now, for 'x' we need to verify, if it is in F_T(S)
          // so we want to check:   F_T^{-1}(x) \in S
          // or alternatively:      G_S^{-1}( F_T^{-1}(x) ) \in G_S^{-1}(S) = T_0

          // we summarize:
          // check if:      G_S^{-1}( F_T^{-1}(x) ) \in T_0

          // S = local_grid_it

          // get geoemetry of the local grid entity:
          const EnGeometryType& geometry_of_S = local_grid_it->geometry();

          // corner of the local element S:
          // ( map the reference corners to the global corners: )
          DomainType corner_0_of_S = geometry_of_S.global(reference_corner_0);
          DomainType corner_1_of_S = geometry_of_S.global(reference_corner_1);
          DomainType corner_2_of_S = geometry_of_S.global(reference_corner_2);

          // Wir bilden die Referenzabbildung G_S : T_0 -> S, gegeben durch
          // G_S( y ) = A_S y + b_S
          // und überprüfen ob G_S^(-1)( F_T^{-1}(x) ) in T_0 liegt

          // value of the matrix A_S (in G(x) = A_S x + a_0) (belonging to the element S)
          double val_A_S[2][2];
          val_A_S[0][0] = corner_1_of_S[0] - corner_0_of_S[0];
          val_A_S[0][1] = corner_2_of_S[0] - corner_0_of_S[0];
          val_A_S[1][0] = corner_1_of_S[1] - corner_0_of_S[1];
          val_A_S[1][1] = corner_2_of_S[1] - corner_0_of_S[1];

          // define 'c := (a_1(1) - a_0(1))·(a_2(2) - a_0(2)) - (a_1(2) - a_0(2))·(a_2(1) - a_0(1))
          double c_S = 1.0 / ( (val_A_S[0][0] * val_A_S[1][1]) - (val_A_S[0][1] * val_A_S[1][0]) );
          // then the inverse A^{-1} is given by:
          // A^{-1}_11 = (1/c) (a_2(2) - a_0(2))     A^{-1}_12 = (1/c) (a_0(1) - a_2(1))
          // A^{-1}_21 = (1/c) (a_0(2) - a_1(2))     A^{-1}_22 = (1/c) (a_1(1) - a_0(1))

          // (A_S^{-1})^T:
          double val_A_S_inverse[2][2];
          val_A_S_inverse[0][0] = c_S * val_A_S[1][1];
          val_A_S_inverse[1][0] = c_S * (-1.0) * val_A_S[1][0];
          val_A_S_inverse[0][1] = c_S * (-1.0) * val_A_S[0][1];
          val_A_S_inverse[1][1] = c_S * val_A_S[0][0];

          // x = F_T( G_S(y) )
          // =>
          // F_T^{-1}(x) = G_S(y) = A_S y + (a_S)_0
          // =>
          // y = A_S^{-1} ( F_T^{-1}(x) - (a_S)_0 )

          // of this 'y', we need to check, if it belongs to T_0

          DomainType G_inverse_of_F_inverse_of_x(0.0);
          for (int k = 0; k < 2; ++k)
            for (int l = 0; l < 2; ++l)
              G_inverse_of_F_inverse_of_x[k] += val_A_S_inverse[k][l] * (F_inverse_of_x[l] - corner_0_of_S[l]);

          // y \in T_0 ?
          if ( ( (G_inverse_of_F_inverse_of_x[0] >= 0.0) && (G_inverse_of_F_inverse_of_x[0] <= 1.0) )
               && ( (G_inverse_of_F_inverse_of_x[1] >= 0.0)
                    && ( G_inverse_of_F_inverse_of_x[1] <= (1.0 - G_inverse_of_F_inverse_of_x[0]) ) ) )
          {
            LocalFunctionType local_corrector_disc_func = corrector_of_coarse_disc_func.localFunction(*local_grid_it);

            // create at quadrature with 3 quadrature points:
            Fem::CachingQuadrature< GridPartType, 0 > local_quad(*local_grid_it, 2);          // 3 points for linear pol in
                                                                                         // 2D
            // for the Lagrange Interpolation on the relevant coarse entity:
            DomainType fine_quad_point[3];
            RangeType local_value_corrector_func[3];

            for (int qp = 0; qp < 3; ++qp)
            {
              fine_quad_point[qp] = geometry_of_S.global( local_quad.point(qp) );

              // evaluate local coarse disc function
              local_corrector_disc_func.evaluate(local_quad[qp], local_value_corrector_func[qp]);
            }

            LinearLagrangeFunction2D< DiscreteFunctionSpaceType >
            corrector_disc_func_entity(fine_quad_point[0],
                                       local_value_corrector_func[0],
                                       fine_quad_point[1],
                                       local_value_corrector_func[1],
                                       fine_quad_point[2],
                                       local_value_corrector_func[2]);

            RangeType corrector_value = 0.0;
            corrector_disc_func_entity.evaluate(x, corrector_value);
            value += corrector_value;

            break;
          }
        }

        break;
      }
    }  // end coarse grid element iteration

    if (error_in_compuation)
    {
      return 0.0;
    }

    return value;
  } // end method


  // for two discrete functions...
  template< typename LocalProblemNumManagerType, int polOrd >
  RangeFieldType error_L2_with_corrector(const DiscreteFunctionType& fine_discrete_function,
                                         const DiscreteFunctionType& coarse_discrete_function /* that can be
                                                                                                 *reconstructed */,
                                         const LocalProblemNumManagerType& lp_num_manager) const {
    const DiscreteFunctionSpaceType& fine_space = fine_discrete_function.space();
    const GridPartType& fineGridPart = fine_space.gridPart();
    const auto& comm = fineGridPart.grid().comm();

    RangeFieldType ret = 0;

    if (dimRange > 1)
    {
      DUNE_THROW(Dune::NotImplemented,"L2Error::norm2: only implemented for dimRange = 1! \n");
    }

    // get the local discrete function space
    const DiscreteFunctionSpaceType& local_space = lp_num_manager.get_local_discrete_function_space();

    // for product:
    // int quadOrd = (polOrd*dim)*(polOrd*dim);
    const int quadOrd = polOrd;

    // iterate over all elements defining the function
    const IteratorType fine_it_end = fine_space.end();
    for (IteratorType fine_it = fine_space.begin(); fine_it != fine_it_end; ++fine_it)
    {
      const EntityType& fine_en = *fine_it;

      const Fem::CachingQuadrature< GridPartType, 0 > quad(fine_en, quadOrd);
      // get local functions on current element
      const LocalFunctionType local_fine_discrete_function = fine_discrete_function.localFunction(fine_en);

      // get geoemetry of entity
      const EnGeometryType& fine_geo = fine_en.geometry();

      const int quadNop = quad.nop();
      for (int qp = 0; qp < quadNop; ++qp)
      {
        const double det
          = fine_geo.integrationElement( quad.point(qp) );

        const DomainType x = fine_geo.global( quad.point(qp) );

        RangeType val_fine_disc_func = 0.0;
        // evaluate local function
        local_fine_discrete_function.evaluate(quad[qp], val_fine_disc_func);

        RangeType val_corrected_coarse_disc_func = 0.0;

        val_corrected_coarse_disc_func
          = evaluate_with_corrector< LocalProblemNumManagerType, 2* DiscreteFunctionSpaceType::polynomialOrder + 2 >(
          coarse_discrete_function,
          x,
          local_space,
          lp_num_manager);

        ret += det * quad.weight(qp)
               * ( (val_fine_disc_func
                    - val_corrected_coarse_disc_func) * (val_fine_disc_func - val_corrected_coarse_disc_func) );
      }     // end qp iteration
    }   // end element iteration

    ret = comm.sum(ret);
    return sqrt(ret);
  } // end method


  // expensive hack:
  template< int polOrd >
  RangeFieldType norm_adaptive_grids(const DiscreteFunctionType& coarse_disc_func,
                                     const DiscreteFunctionType& fine_disc_func, double /*dummy*/ = 0) const {
    // check if the discrete functions have valid dofs:
    if ( !coarse_disc_func.dofsValid() || !fine_disc_func.dofsValid() )
    {
      DSC_LOG_ERROR << "Solution of discrete function invalid." << std::endl;
      return 0.0;
    }
    // get function spaces
    const DiscreteFunctionSpaceType& fine_discreteFunctionSpace = fine_disc_func.space();
    const GridPartType& fine_gridPart = fine_discreteFunctionSpace.gridPart();
    const auto& fine_comm = fine_gridPart.grid().comm();

    // to return the L2 Norm:
    RangeFieldType l2Norm = 0.0;

    // for product:
    // last entity of fine grid:
    const IteratorType fine_end = fine_discreteFunctionSpace.end();
    for (IteratorType fine_it = fine_discreteFunctionSpace.begin(); fine_it != fine_end; ++fine_it)
    {
      // get geoemetry of fine grid entity:
      const EnGeometryType& fine_geo = fine_it->geometry();

      // entity
      const EntityType& fine_entity = *fine_it;
      const Fem::CachingQuadrature< GridPartType, 0 > quadrature(fine_entity, polOrd);
      // integrate
      const int quadratureNop = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint)
      {
        const double det = quadrature.weight(quadraturePoint)
                           * fine_geo.integrationElement( quadrature.point(quadraturePoint) );

        // das Zentrum des Fine-Grid Elements:
        const DomainType global_quad_point = fine_geo.global( quadrature.point(quadraturePoint) );

        RangeType coarse_value;
        coarse_disc_func.evaluate(global_quad_point, coarse_value);

        RangeType fine_value;
        fine_disc_func.evaluate(global_quad_point, fine_value);

        l2Norm += det * pow(coarse_value - fine_value, 2.0);
      }
    }  // end fine grid element iteration

    l2Norm = fine_comm.sum(l2Norm);

    return sqrt(l2Norm);
  } // end method

  // for the following method, the DiscreteFunction 'fine_disc_func' needs to be defined on a gridPart that is a
  // refinment of the gridPart of 'coarse_disc_func'
  template< int polOrd >
  RangeFieldType norm2(const DiscreteFunctionType& coarse_disc_func,
                       const DiscreteFunctionType& fine_disc_func, double /*dummy*/ = 0) const {
    // get function spaces
    const DiscreteFunctionSpaceType& coarse_discreteFunctionSpace = coarse_disc_func.space();
    const DiscreteFunctionSpaceType& fine_discreteFunctionSpace = fine_disc_func.space();
    const GridPartType& fine_gridPart = fine_discreteFunctionSpace.gridPart();
    const GridType& fine_grid = fine_gridPart.grid();
    const auto& fine_comm = fine_gridPart.grid().comm();

    // to return the L2 Norm:
    RangeFieldType l2Norm = 0.0;

    if (dimRange > 1)
    {
      DUNE_THROW(Dune::NotImplemented,"L2Error::norm2: only implemented for dimRange = 1! \n");
    }

    // for product:
    // int quadOrd = (polOrd*2)*(polOrd*2);
    const int quadOrd = polOrd;

    const IdSet& idset_fine_grid = fine_grid.localIdSet();
    int number_of_coarse_elements_that_are_fathers = 0;

    // number of visited coarse grid elements
    int number_of_coarse_grid_elements = 0;

    // last entity of coarse grid:
    const IteratorType coarse_entity_end = coarse_discreteFunctionSpace.end();
    for (IteratorType coarse_it = coarse_discreteFunctionSpace.begin(); coarse_it != coarse_entity_end; ++coarse_it)
    {
      const int coarse_entity_level = coarse_it->level();
      const LocalFunctionType local_coarse_disc_func = coarse_disc_func.localFunction(*coarse_it);

      // create at quadrature with 3 quadrature points:
      const Fem::CachingQuadrature< GridPartType, 0 > coarse_quad(*coarse_it, 2);    // 3 points for linear pol in 2D
      // get geoemetry of coarse entity
      const EnGeometryType& coarse_geo = coarse_it->geometry();
      DomainType coarse_quad_point[3];
      RangeType local_value_coarse_func[3];
      for (int qp = 0; qp < 3; ++qp)
      {
        coarse_quad_point[qp] = coarse_geo.global( coarse_quad.point(qp) );

        // evaluate local coarse disc function
        local_coarse_disc_func.evaluate(coarse_quad[qp], local_value_coarse_func[qp]);
      }

      LinearLagrangeFunction2D< DiscreteFunctionSpaceType >
      coarse_disc_func_entity(coarse_quad_point[0],
                              local_value_coarse_func[0],
                              coarse_quad_point[1],
                              local_value_coarse_func[1],
                              coarse_quad_point[2],
                              local_value_coarse_func[2]);

      const IdType id_coarse_it = idset_fine_grid.id(*coarse_it);

      bool has_father = false;
      // last entity of fine grid:
      IteratorType fine_end = fine_discreteFunctionSpace.end();
      for (IteratorType fine_it = fine_discreteFunctionSpace.begin(); fine_it != fine_end; ++fine_it)
      {
        const int fine_entity_level = fine_it->level();
        const int level_difference = fine_entity_level - coarse_entity_level;

        EntityPointerType fine_father_entity = fine_it;
        for (int lev = 0; lev < level_difference; ++lev)
          fine_father_entity = fine_father_entity->father();

        const IdType id_fine_father_entity = idset_fine_grid.id(*fine_father_entity);

        if (id_fine_father_entity == id_coarse_it)
        {
          has_father = true;
          const Fem::CachingQuadrature< GridPartType, 0 > fine_quad(*fine_it, quadOrd);

          // get local functions on current element
          const LocalFunctionType local_fine_disc_func = fine_disc_func.localFunction(*fine_it);

          // get geoemetry of entity
          const EnGeometryType& fine_geo = fine_it->geometry();

          const int fine_quadNop = fine_quad.nop();
          for (int qp = 0; qp < fine_quadNop; ++qp)
          {
            const double det
              = fine_geo.integrationElement( fine_quad.point(qp) );

            // find best quadrature point in the coarse grid quadrature:
            const DomainType fine_quad_point = fine_geo.global( fine_quad.point(qp) );

            RangeType coarse_value;
            coarse_disc_func_entity.evaluate(fine_quad_point, coarse_value);

            RangeType fine_value = 0.0;
            // evaluate fine local function
            //!TODO local function eval is thrown away
            local_fine_disc_func.evaluate(fine_quad[qp], fine_value);
            fine_value = 1.0;        //! -= coarse_value;

            l2Norm += det * fine_quad.weight(qp) * (fine_value * fine_value);
          }        // end qp iteration
        }
      }

      if (has_father)
      { number_of_coarse_elements_that_are_fathers += 1; } else
      { DSC_LOG_DEBUG << "Coarse element number " << number_of_coarse_grid_elements << "is not a father!!!" << std::endl; }

      number_of_coarse_grid_elements += 1;
    } // end element iteration

    DSC_LOG_INFO  << "Number of visited coarse grid elements: " << number_of_coarse_grid_elements << std::endl
                  << "number_of_coarse_elements_that_are_fathers: " << number_of_coarse_elements_that_are_fathers
                  << std::endl;

    l2Norm = fine_comm.sum(l2Norm);
    return sqrt(l2Norm);
  } // end method
}; // end of class L2Error
#endif
} // end namespace DUNE

#endif // DUNE_MULTISCALE_IMPROVED_L2ERROR_HH
