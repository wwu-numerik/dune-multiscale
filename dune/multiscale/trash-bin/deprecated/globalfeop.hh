#ifndef DUNE_FEOPERATOR_HH
#define DUNE_FEOPERATOR_HH

// - Dune includes
#include <dune/common/fmatrix.hh>
#include <dune/grid/common/referenceelements.hh>

// - local includes
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/localoperator.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>

namespace Dune {
/** @defgroup FEOpInterface FEOpInterface
   *  @ingroup DiscreteOperator
   *
   *
   * @{
   */

/*======================================================================*/
/*!
   *  \class FEOpInterface
   *  \brief FEopInterface is the interface for the definition of a finite
   *         element operator.
   *
   *  This interface is an old version, which is kept for the poisson example,
   *  but a more general finite element operator FEOp is available in
   *  dune/fem/operator/feop.hh , which is not derived from an interface.
   */
/*======================================================================*/

template< class DiscFunctionType, class FEOpImp >
class FEOpInterface
  : public Operator< typename DiscFunctionType::DomainFieldType,
                     typename DiscFunctionType::RangeFieldType, DiscFunctionType, DiscFunctionType >
{
public:
/*======================================================================*/
/*!
   *   getLocalMatrix: Interface method that returns the local matrix of the
   *                   finite element operator on an entity
   *
   *   This method has to be provided by derived classes.
   *
   *   \param the entity, the local matrix size and a reference to an
   *          instance of the local matrix implementation
   */
/*======================================================================*/

  template< class EntityType, class LocalMatrixImp >
  void getLocalMatrix(EntityType& entity, const int matSize, LocalMatrixImp& mat) const {
    return asImp().getLocalMatrix(entity, matSize, mat);
  }

protected:
  // Barton-Nackman
  FEOpImp& asImp() { return static_cast< FEOpImp& >(*this); }

  const FEOpImp& asImp() const { return static_cast< const FEOpImp& >(*this); }
};

/** @} end documentation group */

/*======================================================================*/
/*!
   *  \class FEOp
   *  \brief The FEOp class provides one example of a class satisfying the
   *         FEOpInterface.
   *
   *  The MatrixImp class must provide a storage type for the global
   *  operator matrix and some arithmetics and basic functionality. In
   *  particular a print() method and an apply() method are assumed to exist.
   *  Additionally a constructor with three arguments is required as
   *  explained in newEmptyMatrix. An add(row,col,val) method is assumed
   *  to exist.
   *
   *  different operating modes are possible. In case of ASSEMBLED, the
   *  whole
   *  global operator matrix is allocated and completely precomputed by
   *  corresponding methods. In case of ON_THE_FLY, no complete matrix is
   *  allocated, but the matrix-vector multiplication is performed by
   *  on-the-fly computation of the elementwise local matrices.
   */
/*======================================================================*/

template< class DiscFunctionType, class MatrixImp, class FEOpImp >
class FEOp
  : public FEOpInterface< DiscFunctionType, FEOpImp >
    , public LocalOperatorDefault< DiscFunctionType, DiscFunctionType, typename
                                   DiscFunctionType::RangeFieldType, FEOpImp >
{
public:
  // ! Type of matrix storage class used for global operator matrix
  typedef MatrixImp MatrixType;

  // ! Operation mode: global allocation and multiplication or only
  // ! on-the-fly
  enum OpMode { ON_THE_FLY, ASSEMBLED };

  // ! fix the size of the local matrices to some reasonable extent
  enum { maxnumOfBaseFct = 100 };

/*======================================================================*/
/*!
   *   constructor: Initialization of FEOp
   *
   *   Based on an existing instance of a function space the FEOp is
   *   initialized.
   *   Operation mode must be selected as ASSEMBLED or ON_THE_FLY.
   *
   *   ???? The role of isleaf is unclear ?????
   *
   *   \param an instance of the discrete function space, the operator mode
   *          and a leaf-flag
   *
   *   \return the initialized FEOp
   */
/*======================================================================*/

  FEOp(const typename DiscFunctionType::DiscreteFunctionSpaceType& fuspace,
       OpMode opMode = ASSEMBLED, bool leaf = true)
    : functionSpace_(fuspace)
      , matrix_(0)
      , matrix_assembled_(false)
      , arg_(NULL)
      , dest_(NULL)
      , opMode_(opMode)
      , leaf_(leaf) {}

/*======================================================================*/
/*!
   *   destructor: In case of allocation of global operator matrix
   *               it is deallocated.
   */
/*======================================================================*/

  ~FEOp() {
    if (matrix_) delete matrix_;
  }

public:
/*======================================================================*/
/*!
   *   print: print matrix to standard out
   */
/*======================================================================*/

  void print() const {
    if (!this->matrix_assembled_) this->assemble();
    this->matrix_->print(std::cout);
  }

/*======================================================================*/
/*!
   *   myMatrix: return reference to Matrix for oem solvers.
   *
   *   The assembled matrix is returned. That means, if global matrix is not
   *   yet allocated, a new empty matrix is generated. If current global
   *   matrix is not yet assembled, an assembly is initiated.
   *
   *   ??? What happens in case of ON_THE_FLY ???
   *
   *   \return reference to assembled global matrix
   */
/*======================================================================*/

  MatrixType& systemMatrix() const {
    // assert(matrix_assembled_ == true);
    if (!this->matrix_assembled_)
    {
      if (!this->matrix_)
        this->matrix_ = this->newEmptyMatrix();
      this->assemble();
    }
    return *this->matrix_;
  } // systemMatrix

/*======================================================================*/
/*!
   *   hasPcMatrix: return false
   *
   *   ???  What is this method good for ???
   *
   *   \return false
   */
/*======================================================================*/

  bool hasPcMatrix() const { return false; }

/*======================================================================*/
/*!
   *   pcMatrix: call myMatrix
   *
   *   ???? The use and relevance of this method is unclear ????
   *
   *   \return the result of systemMatrix
   */
/*======================================================================*/

  MatrixType& pcMatrix() const { return systemMatrix(); }

/*======================================================================*/
/*!
   *   operator(): application operator
   *
   *   This method only makes sense in case of ASSEMBLED opmode, as a global
   *   apply() on the matrix is called, i.e. a matrix-vector-multiplication
   *   with the vector arg is performed, the result stored in dest. This is
   *   an implicit requirement on the MatrixImp class!
   *
   *   \param the argument for the matrix-vector multiplication and the
   *          storage for the destination vector
   */
/*======================================================================*/

  virtual void operator()(const DiscFunctionType& arg,
                          DiscFunctionType& dest) const {
    // diese Funktion muss vorhanden, liegt aber auf Eis, das heißt hier passiert nichts für unser Beispiel. Sie
    // könnte genauso gut leer sein.
    assert(this->opMode_ == ASSEMBLED);

    if (!matrix_assembled_)
    {
      assemble();
    }
    matrix_->apply(arg, dest);
  } // ()

public:
/*======================================================================*/
/*!
   *   isLeaf: returns true if Leafiterator should be used, else false is
   *   returned
   *
   *   ??? what is this good for? ???
   *
   *   \return the current isleaf_ value
   */
/*======================================================================*/

  bool isLeaf() {
    return leaf_;
  }

protected:
  // ! the corresponding function_space
  const typename DiscFunctionType::DiscreteFunctionSpaceType & functionSpace_;

  // ! pointer to the representing global matrix
  mutable MatrixType* matrix_;

  // ! flag indicating whether the global matrix is assembled
  mutable bool matrix_assembled_;

  // ! pointers to storage of argument and destination
  const DiscFunctionType* arg_;
  DiscFunctionType* dest_;

protected:
  /*!
     *   newEmptyMatrix: allocation of a new global matrix
     *
     *   The Matrixclass is required to have a constructor with the syntax
     *   MatrixType(nrows, ncols, nonzeros_per_row).
     *
     *   \return a pointer to the newly allocated global matrix.
     */
  MatrixType* newEmptyMatrix() const {
    typedef typename DiscFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridType GridType;

    enum { dimension = GridType::dimension };
    enum { polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder };

    return new MatrixType(this->functionSpace_.size(),
                          this->functionSpace_.size(),
                          15 * (dimension - 1) * polynomialOrder);
  } // newEmptyMatrix

  /*!
     *   assemble: perform grid-walkthrough and assemble global matrix
     *
     *   If the matrix storage is
     *   not allocated, new storage is allocated by newEmptyMatrix.
     *   the begin and end iterators are determined and the assembling
     *   of the global matrix initiated by call of assembleOnGrid and
     *   bndCorrectOnGrid. The assemled flag is set.
     */
  void assemble() const {
    if (this->matrix_ == NULL)
      matrix_ = this->newEmptyMatrix();

    typedef typename DiscFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridType GridType;

    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

    {
      FieldMatrix< double, maxnumOfBaseFct, maxnumOfBaseFct > mat;

      IteratorType it = functionSpace_.begin();
      IteratorType endit = functionSpace_.end();

      for ( ; it != endit; ++it)
        assembleOnGrid(*it, mat);
    }

    // ! start some boundary treatment
    {
      IteratorType it = functionSpace_.begin();
      const IteratorType endit = functionSpace_.end();
      for ( ; it != endit; ++it)
        boundaryCorrectOnGrid(*it);
    }
    // ! end boundary treatment

    matrix_assembled_ = true;
  } // assemble

  /*!
     *   assembleOnGrid: perform grid walkthrough and assemble matrix
     *
     *   For each element, the local element matrix is determined into the
     *   given local matrix storage and distributed into the global matrix.
     *   Distribution is performed by an add(row,col,val) method on the
     *   global matrix class.
     *
     *   ??? why shoudl this one not be private? method assemble() is
     *   sufficient for public call. ???
     *
     *   \param start and end iterator and storage for a local matrix
     */
  template< class EntityType, class LocalMatrixImp >
  void assembleOnGrid(const EntityType& entity,
                      LocalMatrixImp& mat) const {
    typedef typename DiscFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    BaseFunctionSetType;

    const BaseFunctionSetType& baseSet = functionSpace_.baseFunctionSet(entity);
    // the BaseSet is always local. It denotes the set of base functions that are not equal to zero on a current entity.
    // Therefor 'numBaseFunctions()' simply yields the number of elements of such a set (not the global number of base
    // functions)
    // example: 1-D, polynimial order = 1 => numBaseFunctions = 2

    const int numBaseFunctions = baseSet.numBaseFunctions();

    // setup matrix
    getLocalMatrix(entity, numBaseFunctions, mat);

    for (int i = 0; i < numBaseFunctions; ++i)
    {
      const int row = functionSpace_.mapToGlobal(entity, i);
      for (int j = 0; j < numBaseFunctions; ++j)
      {
        const int col = functionSpace_.mapToGlobal(entity, j);
        matrix_->add(row, col, mat[i][j]);
      }
    }
  }   // end method

  /*!
     *   boundaryCorrectOnGrid: treatment of Dirichlet-DOFS
     *
     *   delete rows and columns for dirichlet DOFS, setting diagonal
     *   element to 1. Lagrange Basis is implicitly assumed.
     *
     *   \param start and end iterator
     */
  template< class EntityType >
  void boundaryCorrectOnGrid(const EntityType& entity) const {
    typedef typename DiscFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::LagrangePointSetType
    LagrangePointSetType;

    enum { faceCodim = 1 };
    typedef typename GridPartType::IntersectionIteratorType
    IntersectionIteratorType;
    typedef typename LagrangePointSetType::template Codim< faceCodim >
      ::SubEntityIteratorType
    FaceDofIteratorType;

    const DiscreteFunctionSpaceType& discreteFunctionSpace = functionSpace_;
    const GridPartType& gridPart = discreteFunctionSpace.gridPart();

    IntersectionIteratorType it = gridPart.ibegin(entity);
    const IntersectionIteratorType endit = gridPart.iend(entity);
    for ( ; it != endit; ++it)
    {
      if ( !(*it).boundary() )
        continue;

      const LagrangePointSetType& lagrangePointSet
        = discreteFunctionSpace.lagrangePointSet(entity);

      const int face = (*it).indexInInside();
      FaceDofIteratorType faceIt
        = lagrangePointSet.template beginSubEntity< faceCodim >(face);
      const FaceDofIteratorType faceEndIt
        = lagrangePointSet.template endSubEntity< faceCodim >(face);
      for ( ; faceIt != faceEndIt; ++faceIt)
      {
        const unsigned int dof
          = discreteFunctionSpace.mapToGlobal(entity, *faceIt);
        matrix_->kroneckerKill(dof, dof);
      }
    }
  } // boundaryCorrectOnGrid

private:
  // ! operator mode
  OpMode opMode_;

  // ! true if LeafIterator is used, deprecated because the Iterator now
  // ! comes from the space
  bool leaf_;
};
} // end namespace

#endif // ifndef DUNE_FEOPERATOR_HH
