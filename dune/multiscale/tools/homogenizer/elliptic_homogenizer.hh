/**************************************************************************
  **       Title: L2Error class
  **    $RCSfile$
  **   $Revision: 1723 $$Name$
  **       $Date: 2007-06-20 15:20:54 +0000 (Wed, 20 Jun 2007) $
  **   Copyright: GPL $Author: dune community $
  ** Description: L2 error class, which computes the error between a function
  **              and a discrete function. Extracted class from
  **              Roberts poisson-example
  **
  **************************************************************************/

// das alles klappt momentan nur in 2-D!!! Es laesst sich aber sehr einfach verallgemeinern.
// Implementierung funktioniert unter folgenden Voraussetzungen an den (Diffusions-)Tensor:
// 1. Elliptizitätsbedingung
// 2. A(x,y) = A(y)


#ifndef DUNE_HOMOGENIZER_HH
#define DUNE_HOMOGENIZER_HH

// where the quadratures are defined
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/multiscale/tools/assembler/matrix_assembler/elliptic_fem_matrix_assembler.hh>

#include <dune/fem/gridpart/periodicgridpart.hh>

// for data output:
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>

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

template< class FunctionSpaceImp >
class ZeroFunction
  : public Dune::Fem::Function< FunctionSpaceImp, ZeroFunction< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef ZeroFunction< FunctionSpaceType >       ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

public:
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = 0;
  }
};

template< class FunctionSpaceImp >
class MassWeight
  : public Dune::Fem::Function< FunctionSpaceImp, MassWeight< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef MassWeight< FunctionSpaceType >         ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

protected:
  const RangeFieldType lambda_;

public:
  // Constructor
  inline explicit MassWeight(const RangeFieldType lambda)
    : lambda_(lambda)
      // lambda itself was only valid within ZeroFunction, therefore we needed to save its value in lambda_()
  {}

  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y[0] = lambda_;
  }
};

#if 1
// since we need to evaluate A( x, \cdot ) to solve cellproblems (in comparison to A( \cdot, \frac{\cdot}{\epsilon} )
// for the global problem), we must transform the orginal tensor to be able to use a standard FEM operator for solving
// cell problems (otherwise: calling the method evaluate(i,j,x,y) within the matrixassembler would evaluate A^{\epsilon}
// instead of A(x,\cdot) )
template< class FunctionSpaceImp, class TensorImp >
class TransformTensor
  : public Function< FunctionSpaceImp, TransformTensor< FunctionSpaceImp, TensorImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;
  typedef TensorImp        TensorType;

private:
  typedef TransformTensor< FunctionSpaceType, TensorImp > ThisType;
  typedef Function< FunctionSpaceType, ThisType >         BaseType;

public:
  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

  static const int dimDomain = DomainType::dimension;

protected:
  const TensorType& tensor_;

public:
  // Constructor
  inline explicit TransformTensor(const TensorType& tensor)
    : tensor_(tensor)
  {}

  void diffusiveFlux(const DomainType& y,
                     const JacobianRangeType& direction,
                     JacobianRangeType& flux) const {
    Problem::ModelProblemData model_info;
    
    //! EPSILON BESSER AUS DEM PARAMETER-FILE HOLEN!
    const double epsilon = model_info.getEpsilon();

    DomainType new_y;

    for (int i = 0; i < dimDomain; ++i)
      new_y[i] = epsilon * y[i];

    tensor_.diffusiveFlux(new_y, direction, flux);
  } // diffusiveFlux

  inline void evaluate(const int i, const int j,
                       const DomainType& x,
                       const TimeType& time,
                       RangeType& y) const {
    std::cout << "Do not use this evaluate()-method !!" << std::endl;
    abort();
    evaluate(i, j, x, y);
  }

  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    std::cout << "Do not use this evaluate()-method !!" << std::endl;
    abort();
    y = 0;
  }

  inline void evaluate(const DomainType& x,
                       const TimeType& time,
                       RangeType& y) const {
    std::cout << "Do not use this evaluate()-method !!" << std::endl;
    abort();
    y = 0;
  }
};
#endif // if 1

// the following class is comparable to a SecondSource-Class (some kind of -div G )
template< class FunctionSpaceImp, class TensorImp >
class CellSource
  : public Dune::Fem::Function< FunctionSpaceImp, CellSource< FunctionSpaceImp, TensorImp > >
{
public:
  typedef TensorImp TensorType;

  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef CellSource< FunctionSpaceType, TensorImp > ThisType;
  typedef Function< FunctionSpaceType, ThisType >    BaseType;

public:
  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

public:
  const FunctionSpaceType& functionSpace_;
  const TensorType& tensor_;
  const int& j_;

  inline explicit CellSource(const FunctionSpaceType& functionSpace, const TensorType& tensor, const int& j)
    : functionSpace_(functionSpace)
      , tensor_(tensor)
      , j_(j) // we solve the j'th cell problem
  {}

  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y[0] = 0;
  }

  inline void evaluate(const int i, const DomainType& y,
                       RangeType& z) const {
    JacobianRangeType direction;
    JacobianRangeType flux;
    
    for (int i_ = 0; i_ < dimDomain; ++i_)
    {
       if (j_ == i_)
       {
	 direction[0][i_] = 1.0;
       }
       else
       {
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
  : public Dune::Fem::Function< FunctionSpaceImp, DefaultDummyAdvection< FunctionSpaceImp > >
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

private:
  typedef DefaultDummyAdvection< FunctionSpaceType > ThisType;
  typedef Function< FunctionSpaceType, ThisType >    BaseType;

public:
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType  RangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

  const FunctionSpaceType* functionSpace_;

protected:
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

  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       const DomainType& y,
                       RangeType& z) const {
    z = 0;
  }

  inline void evaluate(const int i,
                       const DomainType& x,
                       const DomainType& y,
                       RangeType& z) const {
    z = 0;
  }

  inline void evaluate(const int i,
                       const int j,
                       const DomainType& x,
                       RangeType& y) const {
    y = 0;
  }

  inline void evaluate(const int i,
                       const DomainType& x,
                       RangeType& y) const {
    y = 0;
  }

  inline void evaluate(const int i,
                       const DomainType& x,
                       const TimeType& t,
                       RangeType& y) const {
    y = 0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       RangeType& y) const {
    y = 0;
  }

  // dummy implementation
  inline void evaluate(const DomainType& x,
                       const TimeType time,
                       RangeType& y) const {
    y = 0;
  }
};

template< class GridImp, class TensorImp >
class Homogenizer
{
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
  EntityPointerType; // !Brauchen wie das? Loeschen? loeschen?

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
  std::string& filename_;

public:
  Homogenizer(std::string& filename)
    : filename_(filename)
  {}

private:
  double getEntry(TransformTensorType& tensor,
                  PeriodicDiscreteFunctionSpaceType& periodicDiscreteFunctionSpace,
                  PeriodicDiscreteFunctionType& w_i,
                  PeriodicDiscreteFunctionType& w_j,
                  const int& i,
                  const int& j) {
    double a_ij_hom = 0;

    DomainType dummy(0.0);

    IteratorType endit = periodicDiscreteFunctionSpace.end();

    for (IteratorType it = periodicDiscreteFunctionSpace.begin(); it != endit; ++it)
    {
      // entity
      const EntityType& entity = *it;

      PeriodicLocalFunctionType localW_i = w_i.localFunction(entity);
      PeriodicLocalFunctionType localW_j = w_j.localFunction(entity);

      // create quadrature for given geometry type
      QuadratureType quadrature(entity, 2);

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
          // redundant (just for tests):
          #if 0
          direction_j[0][k] += grad_w_j[0][k];
          #endif
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

  FieldMatrix< RangeType, dimension, dimension > getHomTensor(TensorType& tensor) {
    FieldMatrix< RangeType, dimension, dimension > a_hom;

    // to solve cell problems, we always need to use a perforated unit cube as domain:
    GridPtr< GridType > periodicgridptr(filename_);

    periodicgridptr->globalRefine(10);

    PeriodicGridPartType periodicGridPart(*periodicgridptr);

    PeriodicDiscreteFunctionSpaceType periodicDiscreteFunctionSpace(periodicGridPart);

    // to avoid confusions:
    DummySpaceType dummySpace(periodicGridPart);
    // (sometimes periodicDiscreteFunctionSpace is only a dummy)

    // ! define the type of the corresponding solutions ( discrete functions of the type 'DiscreteFunctionType'):

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

    // define mass (just for cell problems \lambda w - \div A \nabla w = rhs)
    MassWeightType mass(lambda);

    TransformTensorType tensor_transformed(tensor);

    // if we have some additional source term (-div G), define:
    CellSourceType G_0(periodicDiscreteFunctionSpace, tensor_transformed, 0);   // 0'th cell problem
    CellSourceType G_1(periodicDiscreteFunctionSpace, tensor_transformed, 1);   // 1'th cell problem
    // - div ( A \nabla u^{\epsilon} ) = f - div G

    // quite a dummy. It's always f = 0
    const ZeroFunctionType zero;

    // ! build the left hand side (lhs) of the problem

    EllipticOperatorType discrete_cell_elliptic_op(periodicDiscreteFunctionSpace, tensor_transformed, mass);

    FEMMatrix lhsMatrix("Cell Problem Stiffness Matrix", periodicDiscreteFunctionSpace, periodicDiscreteFunctionSpace);
    discrete_cell_elliptic_op.assemble_matrix(lhsMatrix, false /*no boundary treatment*/);

    // ! build the right hand side (rhs) of the problem

    // the same right hand side for HM and FEM methods:
    RightHandSideAssembler< PeriodicDiscreteFunctionType > // !, TransformTensorType, CellSourceType >
    cell_0_assembler;    // !( transTensor, G_0 );

    RightHandSideAssembler< PeriodicDiscreteFunctionType > // !, TransformTensorType, CellSourceType >
    cell_1_assembler;    // !( transTensor, G_1 );

    // Alternativly it is possible to call the RightHandSideAssembler with a second source Term '- div G':
    // RightHandSideAssembler< DiscreteFunctionType > rhsassembler( tensor , G );
    cell_0_assembler.template assemble< 2* PeriodicDiscreteFunctionSpaceType::polynomialOrder >(zero, G_0, rhs_0);       //
                                                                                                                         //
                                                                                                                         //
                                                                                                                         //!
                                                                                                                         //
                                                                                                                         //
                                                                                                                         //G_0
                                                                                                                         //
                                                                                                                         //
                                                                                                                         //loeschen!

    cell_1_assembler.template assemble< 2* PeriodicDiscreteFunctionSpaceType::polynomialOrder >(zero, G_1, rhs_1);       //
                                                                                                                         //
                                                                                                                         //
                                                                                                                         //!
                                                                                                                         //
                                                                                                                         //
                                                                                                                         //G_1
                                                                                                                         //
                                                                                                                         //
                                                                                                                         //loeschen!

    // solve the linear systems (with Bi-CG):

    InverseFEMMatrix fembiCG(lhsMatrix, 1e-8, 1e-8, 20000, VERBOSE);

    fembiCG(rhs_0, cellSolution_0);
    fembiCG(rhs_1, cellSolution_1);

    a_hom[0][0] = getEntry(tensor_transformed, periodicDiscreteFunctionSpace, cellSolution_0, cellSolution_0, 0, 0);
    a_hom[1][1] = getEntry(tensor_transformed, periodicDiscreteFunctionSpace, cellSolution_1, cellSolution_1, 1, 1);
    a_hom[0][1] = getEntry(tensor_transformed, periodicDiscreteFunctionSpace, cellSolution_0, cellSolution_1, 0, 1);
    a_hom[1][0] = a_hom[0][1];

    std::cout << "A_homogenized[0][0] = " << a_hom[0][0] << std::endl;
    std::cout << "A_homogenized[0][1] = " << a_hom[0][1] << std::endl;
    std::cout << "A_homogenized[1][0] = " << a_hom[1][0] << std::endl;
    std::cout << "A_homogenized[1][1] = " << a_hom[1][1] << std::endl;

    // in case you want to save the solutions of the two cell problems:
    #if 0
    typedef Tuple< PeriodicDiscreteFunctionType* > IOTupleType;
    typedef DataOutput< GridType, IOTupleType >    DataOutputType;

    // general output parameters
    CellDataOutputParameters outputparam;

    // sequence stamp
    std::stringstream outstring;

    // ------- cell problem 0 -------------

    // create and initialize output class
    IOTupleType cellproblem_0_tuple(&cellSolution_0);

    outputparam.set_prefix("cellSolution_0_");
    DataOutputType cellSolution_0_dataoutput(periodicGridPart.grid(), cellproblem_0_tuple, outputparam);

    // write data
    outstring << "cellSolution_0_";
    cellSolution_0_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
    // clear the std::stringstream:
    outstring.str( std::string() );

    // ------- cell problem 1 -------------

    // create and initialize output class
    IOTupleType cellproblem_1_tuple(&cellSolution_1);

    outputparam.set_prefix("cellSolution_1_");
    DataOutputType cellSolution_1_dataoutput(periodicGridPart.grid(), cellproblem_1_tuple, outputparam);

    // write data
    outstring << "cellSolution_1_";
    cellSolution_1_dataoutput.writeData( 1.0 /*dummy*/, outstring.str() );
    // clear the std::stringstream:
    outstring.str( std::string() );

    // ------------------------------------
    #endif // if 0

    return a_hom;
  } // getHomTensor
}; // end of class
} // end namespace
#endif // ifndef DUNE_HOMOGENIZER_HH
