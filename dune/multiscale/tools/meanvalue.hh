// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_MS_MEANVALUE_HH
#define DUNE_MS_MEANVALUE_HH


#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/misc/l2error.hh>
#include <dune/stuff/fem/functions/integrals.hh>

#include "misc/linear-lagrange-interpolation.hh"

namespace Dune {
/**  \class Meanvalue
   *  \brief The Meanvalue class provides a method to calculate the meanvalue of a discrete function
   *
   *  Actually only the meanvalue of discrete functions on the unit-cube can be calculated.
   *  If you want more, divide the return value of getMeanvalue() by the size of the domain.
   *  Since it's the unit cube for our purpose, the domain-size is 1 and therefore unimportant
   *  for us.
   *
 * Usage of the meanvalue class:
  * Example:
   *
   *  Meanvalue< DiscreteFunctionType > mymean;
   *  DiscreteFunctionSpaceType :: RangeType theMeanValue;
   *  theMeanValue = mymean.getMeanvalue( a_discreteFunction );
   *  std :: cout << "Meanvalue of the numerical solution: " << theMeanValue << std :: endl;
   *
   * Shift the discrete function to meanvalue zero:
   *
   * mymean.adapt<DiscreteFunctionType>( a_discreteFunction , theMeanValue );
   *
   * \TODO use stuff/fem/funtions/integrals.hh instead
   **/
template <class DiscreteFunctionType>
class Meanvalue {
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::RangeType RangeType;
  typedef typename DiscreteFunctionType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename DiscreteFunctionSpaceType::GridType GridType;
  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef DomainFieldType TimeType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename GridType::template Codim<0>::Entity EntityType;
  typedef typename GridType::template Codim<0>::Geometry EnGeometryType;
  typedef typename EntityType::ctype coordType;

  static const int dimension = GridType::dimension;
  static const int spacePolOrd = DiscreteFunctionSpaceType::CommonTraits::polynomial_order;

  struct FunctorBase {
    virtual void evaluate(const DomainType& global_point, RangeType& y) const = 0;
  };

  RangeType getMeanvalue_common(const DiscreteFunctionSpaceType& space, const FunctorBase& function) const {
    RangeType y(0.0);
    RangeType theMeanValue(0.0);

    for (const auto& entity : space) {
      const auto quadrature = DSFe::make_quadrature(entity, space);
      const auto& geo = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint) {
        const double det =
            quadrature.weight(quadraturePoint) * geo.integrationElement(quadrature.point(quadraturePoint));

        y = function(geo.global(quadrature.point(quadraturePoint)));

        theMeanValue += det * y;
      }
    }
  }

public:
  RangeType getMeanvalue(const DiscreteFunctionType& discFunc) const {
    const DiscreteFunctionSpaceType& space = discFunc.space();
    RangeType y(0.0); // return value
    RangeType theMeanValue(0.0);

    for (const auto& entity : space) {
      // create quadrature for given geometry type
      const auto quadrature = DSFe::make_quadrature(entity, space);
      const auto localfunc = discFunc.localFunction(entity);
      const auto& geo = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.size();
      for (int quadraturePoint = 0; quadraturePoint < quadratureNop; ++quadraturePoint) {
        const double det =
            quadrature.weight(quadraturePoint) * geo.integrationElement(quadrature.point(quadraturePoint));

        localfunc.evaluate(quadrature, quadraturePoint, y);

        theMeanValue += det * y;
      }
    }

    const auto& comm = space.gridPart().grid().comm();
    theMeanValue = comm.sum(theMeanValue);

    return theMeanValue;
  } // end of method

  template <class FunctionType>
  RangeType getMeanvalue(const DiscreteFunctionSpaceType& space, const FunctionType& function) const {
    struct Functor : public FunctorBase {
      const FunctionType& function;
      Functor(const FunctionType& f) : function(f) {}
      virtual RangeType operator()(const DomainType& global_point, RangeType& y) const {
        function.evaluate(global_point, y);
      }
    } f{function};
    return getMeanvalue_common(space, f);
  } // end of method

  // the case the function is a vector (for instance advection)
  template <class FunctionType>
  RangeType getMeanvalue(const DiscreteFunctionSpaceType& space, const FunctionType& function,
                         const int& i /*in case there are several components*/) const {
    struct Functor : public FunctorBase {
      const FunctionType& function;
      const int i;
      Functor(const FunctionType& f, const int _i)
        : function(f)
        , i(_i) {}
      virtual RangeType operator()(const DomainType& global_point, RangeType& y) const {
        function.evaluate(i, global_point, y);
      }
    } f{function, i};
    return getMeanvalue_common(space, f);
  } // end of method

  // the case the function is a time-dependent vector (for instance advection). The Time t is fixed.
  template <class FunctionType>
  RangeType getMeanvalue(const DiscreteFunctionSpaceType& space, const FunctionType& function, const TimeType& t,
                         const int& i /*in case there are several components*/) const {
    struct Functor : public FunctorBase {
      const FunctionType& function;
      const int i;
      const TimeType t;
      Functor(const FunctionType& f, const int _i, const TimeType _t)
        : function(f)
        , i(_i)
        , t(_t) {}
      virtual RangeType operator()(const DomainType& global_point, RangeType& y) const {
        function.evaluate(i, global_point, t, y);
      }
    } f{function, i, t};
    return getMeanvalue_common(space, f);
  } // end of method

  // the case the function is a matrix (for instance diffusion)
  template <class FunctionType>
  RangeType getMeanvalue(const DiscreteFunctionSpaceType& space, const FunctionType& function, const int& i,
                         const int& j) const {
    struct Functor : public FunctorBase {
      const FunctionType& function;
      const int i;
      const int j;
      virtual RangeType operator()(const DomainType& global_point, RangeType& y) const {
        function.evaluate(i, j, global_point, y);
      }
    } f{function, i, j};
    return getMeanvalue_common(space, f);
  } // end of method

  // the case the function is a matrix (for instance diffusion) with Time t
  template <class FunctionType>
  RangeType getMeanvalue(const DiscreteFunctionSpaceType& space, const FunctionType& function, const TimeType& t,
                         const int& i, const int& j) const {
    struct Functor : public FunctorBase {
      const FunctionType& function;
      const int i;
      const int j;
      const TimeType t;
      Functor(const FunctionType& f, const int _i, const int _j, const TimeType _t)
        : function(f)
        , i(_i)
        , j(_j)
        , t(_t) {}
      virtual RangeType operator()(const DomainType& global_point, RangeType& y) const {
        function.evaluate(i, j, global_point, t, y);
      }
    } f{function, i, j, t};
    return getMeanvalue_common(space, f);
  } // end of method

  //! Subdraktion des Mittelwertes von der DiscreteFunction
  template <class FunctionType>
  static void adapt(FunctionType& discreteFunction, RangeType& meanvalue) {
    typedef typename FunctionType::DofIteratorType DofIteratorType;

    const DofIteratorType end = discreteFunction.dend();
    for (DofIteratorType it = discreteFunction.dbegin(); it != end; ++it)
      *it -= meanvalue[0];
    // Dof Iterator verwenden um Verschiebung um Mittelwert

  } // adapt
};  // end of class Meanvalue

} // end namespace DUNE

#endif // ifndef DUNE_MS_MEANVALUE_HH
