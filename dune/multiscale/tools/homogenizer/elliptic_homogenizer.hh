// das alles klappt momentan nur in 2-D!!! Es laesst sich aber sehr einfach verallgemeinern.
// Implementierung funktioniert unter folgenden Voraussetzungen an den (Diffusions-)Tensor:
// 1. Elliptizitätsbedingung
// 2. A(x,y) = A(y)

#ifndef DUNE_HOMOGENIZER_HH
#define DUNE_HOMOGENIZER_HH

#include <dune/multiscale/tools/assembler/matrix_assembler/elliptic_fem_matrix_assembler.hh>
#include <dune/multiscale/grids/periodicgridpart/periodicgridpart.hh>

// for data output:
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/stuff/fem/functions/analytical.hh>
#include <dune/multiscale/tools/assembler/righthandside_assembler.hh>

namespace Dune {
// define output traits
struct CellDataOutputParameters
  : public DataOutputParameters
{
public:
  std::string my_prefix_;

  // path where the data is stored
  std::string path() const {
    return "data_output_hmm";
  }

  void set_prefix(std::string my_prefix) {
    my_prefix_ = my_prefix;
    // std :: cout << "Set prefix. my_prefix_ = " << my_prefix_ << std :: endl;
  }

  // base of file name for data file
  std::string prefix() const {
    if (my_prefix_ == "")
      return "solutions_cell_problems";
    else
      return my_prefix_;
  }

  // format of output:
  int outputformat() const {
    // return 0; // GRAPE (lossless format)
    return 1; // VTK
    // return 2; // VTK vertex data
    // return 3; // gnuplot
  }
};


NULLFUNCTION(ZeroFunction)

template< class FunctionSpaceImp >
class MassWeight
  : public Dune::Stuff::Fem::ConstantFunction< FunctionSpaceImp >
{
public:
  // Constructor
  inline explicit MassWeight(const typename FunctionSpaceImp::RangeType& lambda)
    : Dune::Stuff::Fem::ConstantFunction< FunctionSpaceImp >(lambda)
  {}

};

/** since we need to evaluate A( x, \cdot ) to solve cellproblems (in comparison to A( \cdot, \frac{\cdot}{\epsilon} )
 * for the global problem), we must transform the orginal tensor to be able to use a standard FEM operator for solving
 * cell problems (otherwise: calling the method evaluate(i,j,x,y) within the matrixassembler would evaluate A^{\epsilon}
 * instead of A(x,\cdot) )
 **/
template< class FunctionSpaceImp, class TensorImp >
class TransformTensor
  : public Fem::Function< FunctionSpaceImp, TransformTensor< FunctionSpaceImp, TensorImp > >
{
private:
  typedef FunctionSpaceImp FunctionSpaceType;
  typedef TensorImp        TensorType;


  typedef TransformTensor< FunctionSpaceType, TensorImp > ThisType;
  typedef Fem::Function< FunctionSpaceType, ThisType >         BaseType;

  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

  static const int dimDomain = DomainType::dimension;

  const TensorType& tensor_;

public:
  // Constructor
  inline explicit TransformTensor(const TensorType& tensor)
    : tensor_(tensor)
  {}

  void diffusiveFlux(const DomainType& y,
                     const JacobianRangeType& direction,
                     JacobianRangeType& flux) const {
    const double epsilon = DSC_CONFIG_GET("problem.epsilon", 1.0f);

    DomainType new_y;

    for (int i = 0; i < dimDomain; ++i)
      new_y[i] = epsilon * y[i];

    tensor_.diffusiveFlux(new_y, direction, flux);
  } // diffusiveFlux

  inline void evaluate(const int /*i*/, const int /*j*/,
                       const DomainType& /*x*/,
                       const TimeType& /*time*/,
                       RangeType& /*y*/) const {
    DUNE_THROW(Dune::InvalidStateException,"Do not use this evaluate()-method !!");
  }

  inline void evaluate(const DomainType& /*x*/,
                       RangeType& /*y*/) const {
    DUNE_THROW(Dune::InvalidStateException,"Do not use this evaluate()-method !!");
  }

  inline void evaluate(const DomainType& /*x*/,
                       const TimeType& /*time*/,
                       RangeType& /*y*/) const {
    DUNE_THROW(Dune::InvalidStateException,"Do not use this evaluate()-method !!");
  }
};

// the following class is comparable to a SecondSource-Class (some kind of -div G )
template< class FunctionSpaceImp, class TensorImp >
class CellSource
  : public Fem::Function< FunctionSpaceImp, CellSource< FunctionSpaceImp, TensorImp > >
{
private:
  typedef TensorImp TensorType;

  typedef FunctionSpaceImp FunctionSpaceType;

  typedef CellSource< FunctionSpaceType, TensorImp > ThisType;
  typedef Fem::Function< FunctionSpaceType, ThisType >    BaseType;

  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;


  const FunctionSpaceType& functionSpace_;
  const TensorType& tensor_;
  const int& j_;

public:
  inline explicit CellSource(const FunctionSpaceType& functionSpace, const TensorType& tensor, const int& j)
    : functionSpace_(functionSpace)
      , tensor_(tensor)
      , j_(j) // we solve the j'th cell problem
  {}

  inline void evaluate(const DomainType& /*x*/,
                       RangeType& y) const {
    y[0] = 0;
  }

  inline void evaluate(const int i, const DomainType& y,
                       RangeType& z) const {
    JacobianRangeType direction;
    JacobianRangeType flux;

    for (int i_ = 0; i_ < DomainType::dimension; ++i_)
    {
      if (j_ == i_)
      {
        direction[0][i_] = 1.0;
      } else {
        direction[0][i_] = 0.0;
      }
    }

    tensor_.diffusiveFlux(y, direction, flux);

    // tensor_.evaluate( i, j_, y, z);
    z = -flux[0][i];
  } // evaluate
};

template< class FunctionSpaceImp >
class DefaultDummyAdvection
  : public Fem::Function< FunctionSpaceImp, DefaultDummyAdvection< FunctionSpaceImp > >
{
private:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef DefaultDummyAdvection< FunctionSpaceType > ThisType;
  typedef Fem::Function< FunctionSpaceType, ThisType >    BaseType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

  const FunctionSpaceType* functionSpace_;
  const double epsilon_;

public:
  // Constructor for Dummy
  inline explicit DefaultDummyAdvection(const FunctionSpaceType& functionSpace)
    : functionSpace_(&functionSpace)
      , epsilon_(0)
  {}

  inline explicit DefaultDummyAdvection(const FunctionSpaceType& functionSpace, const double& epsilon)
    : functionSpace_(&functionSpace)
      , epsilon_(epsilon)
  {}

  inline void evaluate(const int /*i*/,
                       const int /*j*/,
                       const DomainType& /*x*/,
                       const DomainType& /*y*/,
                       RangeType& z) const {
    z = 0;
  }

  inline void evaluate(const int /*i*/,
                       const DomainType& /*x*/,
                       const DomainType& /*y*/,
                       RangeType& z) const {
    z = 0;
  }

  inline void evaluate(const int /*i*/,
                       const int /*j*/,
                       const DomainType& /*x*/,
                       RangeType& y) const {
    y = 0;
  }

  inline void evaluate(const int /*i*/,
                       const DomainType& /*x*/,
                       RangeType& y) const {
    y = 0;
  }

  inline void evaluate(const int /*i*/,
                       const DomainType& /*x*/,
                       const TimeType& /*t*/,
                       RangeType& y) const {
    y = 0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& /*x*/,
                       RangeType& y) const {
    y = 0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& /*x*/,
                       const TimeType /*time*/,
                       RangeType& y) const {
    y = 0;
  }
};

template< class GridImp, class TensorImp >
class Homogenizer
{
private:
  typedef GridImp GridType;
  enum { dimension = GridType::dimension };

  typedef TensorImp TensorType;

  typedef PeriodicLeafGridPart< GridType > PeriodicGridPartType;

  typedef FunctionSpace< double, double, dimension, 1 > FunctionSpaceType;

  typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, PeriodicGridPartType, 1 >
  PeriodicDiscreteFunctionSpaceType;

  // to avoid confusion:
  typedef PeriodicDiscreteFunctionSpaceType DummySpaceType;
  // (sometimes PeriodicDiscreteFunctionSpaceType is only a dummy)

  typedef AdaptiveDiscreteFunction< PeriodicDiscreteFunctionSpaceType > PeriodicDiscreteFunctionType;

  // to avoid confusion:
  typedef PeriodicDiscreteFunctionType DummyType;
  // (sometimes PeriodicDiscreteFunctionType is only a dummy)

  typedef MassWeight< FunctionSpaceType > MassWeightType;

  typedef ZeroFunction< FunctionSpaceType > ZeroFunctionType;

  typedef DefaultDummyAdvection< FunctionSpaceType > DefaultDummyAdvectionType;

  typedef TransformTensor< FunctionSpaceType, TensorType > TransformTensorType;

  typedef CellSource< FunctionSpaceType, TransformTensorType > CellSourceType;

  typedef typename PeriodicDiscreteFunctionSpaceType::JacobianRangeType
  PeriodicJacobianRangeType;

  typedef typename PeriodicDiscreteFunctionType::LocalFunctionType PeriodicLocalFunctionType;

  typedef CachingQuadrature< PeriodicGridPartType, 0 > QuadratureType;

  typedef typename GridType::template Codim< 0 >::Geometry GeometryType;

  typedef typename PeriodicDiscreteFunctionSpaceType::IteratorType IteratorType;

  typedef typename GridType::template Codim< 0 >::Entity EntityType;

  typedef typename GridType::template Codim< 0 >::EntityPointer
  EntityPointerType; //!Brauchen wie das? Loeschen? loeschen?

  struct MatrixTraits
  {
    typedef PeriodicDiscreteFunctionSpaceType                          RowSpaceType;
    typedef PeriodicDiscreteFunctionSpaceType                          ColumnSpaceType;
    typedef LagrangeMatrixSetup< false >                               StencilType;
    typedef ParallelScalarProduct< PeriodicDiscreteFunctionSpaceType > ParallelScalarProductType;

    template< class M >
    struct Adapter
    {
      typedef LagrangeParallelMatrixAdapter< M > MatrixAdapterType;
    };
  };

  typedef SparseRowMatrixOperator< PeriodicDiscreteFunctionType, PeriodicDiscreteFunctionType, MatrixTraits > FEMMatrix;

  typedef OEMBICGSTABOp< PeriodicDiscreteFunctionType, FEMMatrix > InverseFEMMatrix;

  // discrete elliptic operator (corresponds with FEM Matrix)
  typedef DiscreteEllipticOperator< PeriodicDiscreteFunctionType, TransformTensorType,
                                    MassWeightType > EllipticOperatorType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  enum { spacePolOrd = PeriodicDiscreteFunctionSpaceType::polynomialOrder };

  // dgf file that describes the perforated domain
  const std::string& filename_;

public:
  Homogenizer(const std::string& filename)
    : filename_(filename)
  {}

  typedef FieldMatrix< RangeType, dimension, dimension > HomTensorType;
private:
  double getEntry(const TransformTensorType& tensor,
                  const PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                  const PeriodicDiscreteFunctionType& w_i,
                  const PeriodicDiscreteFunctionType& w_j,
                  const int& i,
                  const int& j) const {
    double a_ij_hom = 0;

    DomainType dummy(0.0);

    const IteratorType endit = periodicDiscreteFunctionSpace.end();

    for (IteratorType it = periodicDiscreteFunctionSpace.begin(); it != endit; ++it)
    {
      // entity
      const EntityType& entity = *it;

      PeriodicLocalFunctionType localW_i = w_i.localFunction(entity);
      PeriodicLocalFunctionType localW_j = w_j.localFunction(entity);

      // create quadrature for given geometry type
      const QuadratureType quadrature(entity, 2);

      // get geoemetry of entity
      const GeometryType& geometry = entity.geometry();

      // integrate
      const int quadratureNop = quadrature.nop();
      for (int localQuadPoint = 0; localQuadPoint < quadratureNop; ++localQuadPoint)
      {
        RangeType localIntegral = 0;

        PeriodicJacobianRangeType grad_w_i;
        localW_i.jacobian(quadrature[localQuadPoint], grad_w_i);

        PeriodicJacobianRangeType grad_w_j;
        localW_j.jacobian(quadrature[localQuadPoint], grad_w_j);

        // local (barycentric) coordinates (with respect to cell grid entity)
        const typename QuadratureType::CoordinateType& local_point = quadrature.point(localQuadPoint);

        // global point in the unit cell Y
        const DomainType global_point_in_Y = geometry.global(local_point);

        PeriodicJacobianRangeType direction_i;
        for (int k = 0; k < dimension; ++k)
        {
          direction_i[0][k] = grad_w_i[0][k];
          if (k == i)
          { direction_i[0][k] += 1.0; }
        }

        PeriodicJacobianRangeType direction_j;
        for (int k = 0; k < dimension; ++k)
        {
          direction_j[0][k] = 0.0;
          if (k == j)
          { direction_j[0][k] += 1.0; }
        }

        PeriodicJacobianRangeType flux_i;
        tensor.diffusiveFlux(global_point_in_Y, direction_i, flux_i);

        localIntegral = (flux_i[0] * direction_j[0]);

        const double entityVolume = quadrature.weight(localQuadPoint)
                                    * geometry.integrationElement( quadrature.point(localQuadPoint) );

        a_ij_hom += entityVolume * localIntegral;
      }
    }

    return a_ij_hom;
  } // end of method

public:
  HomTensorType getHomTensor(const TensorType& tensor) const {
    HomTensorType a_hom;

    // to solve cell problems, we always need to use a perforated unit cube as domain:
    GridPtr< GridType > periodicgridptr(filename_);

    periodicgridptr->globalRefine(10);

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
    const CellSourceType G_0(periodicDiscreteFunctionSpace, tensor_transformed, 0);   // 0'th cell problem
    const CellSourceType G_1(periodicDiscreteFunctionSpace, tensor_transformed, 1);   // 1'th cell problem
    // - div ( A \nabla u^{\epsilon} ) = f - div G

    // quite a dummy. It's always f = 0
    const ZeroFunctionType zero;

    //! build the left hand side (lhs) of the problem

    // define mass (just for cell problems \lambda w - \div A \nabla w = rhs)
    const EllipticOperatorType discrete_cell_elliptic_op(periodicDiscreteFunctionSpace,
                                                         tensor_transformed,
                                                         new MassWeightType(lambda));

    FEMMatrix lhsMatrix("Cell Problem Stiffness Matrix", periodicDiscreteFunctionSpace, periodicDiscreteFunctionSpace);
    discrete_cell_elliptic_op.assemble_matrix(lhsMatrix, false /*no boundary treatment*/);

    //! build the right hand side (rhs) of the problem

    // the same right hand side for HM and FEM methods:
    typedef RightHandSideAssembler< PeriodicDiscreteFunctionType > RhsAssembler;

    // Alternativly it is possible to call the RightHandSideAssembler with a second source Term '- div G':
    // RightHandSideAssembler< DiscreteFunctionType > rhsassembler( tensor , G );
    RhsAssembler::template assemble< 2* PeriodicDiscreteFunctionSpaceType::polynomialOrder >(zero, G_0, rhs_0);
    RhsAssembler::template assemble< 2* PeriodicDiscreteFunctionSpaceType::polynomialOrder >(zero, G_1, rhs_1);
    // solve the linear systems (with Bi-CG):

    const InverseFEMMatrix fembiCG(lhsMatrix, 1e-8, 1e-8, 20000, DSC_CONFIG_GET("global.cgsolver_verbose", false));

    fembiCG(rhs_0, cellSolution_0);
    fembiCG(rhs_1, cellSolution_1);

    a_hom[0][0] = getEntry(tensor_transformed, periodicDiscreteFunctionSpace, cellSolution_0, cellSolution_0, 0, 0);
    a_hom[1][1] = getEntry(tensor_transformed, periodicDiscreteFunctionSpace, cellSolution_1, cellSolution_1, 1, 1);
    a_hom[0][1] = getEntry(tensor_transformed, periodicDiscreteFunctionSpace, cellSolution_0, cellSolution_1, 0, 1);
    a_hom[1][0] = a_hom[0][1];

    DSC_LOG_DEBUG << "A_homogenized[0][0] = " << a_hom[0][0] << std::endl
                  << "A_homogenized[0][1] = " << a_hom[0][1] << std::endl
                  << "A_homogenized[1][0] = " << a_hom[1][0] << std::endl
                  << "A_homogenized[1][1] = " << a_hom[1][1] << std::endl;

    return a_hom;
  } // getHomTensor
}; // end of class
} // end namespace
#endif // ifndef DUNE_HOMOGENIZER_HH
