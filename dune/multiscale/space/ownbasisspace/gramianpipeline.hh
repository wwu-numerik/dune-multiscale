#ifndef  DUNE_RB_GRAMIANPIPELINE_HH__
#define  DUNE_RB_GRAMIANPIPELINE_HH__

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include "../discfunclist/wrapper.hh"

namespace Dune {

  namespace Multiscale {

    /** @brief GramianPipeline efficiently computes gramian matrices of a given
     * discrete function list, "gramian vectors" between the functions of a
     * discrete function list and a given discrete function. The discrete functions
     * from the list can also be applied to discrete operators first.
     *
     * Assuming you have
     *   - a discrete function space @f${\cal W}_h@f$,
     *   - a list of discrete function snapshots @f$\Phi:=\{\varphi^k\}_{k=1}^K@f$,
     *     with @f$\varphi^k \in {\cal W}_h@f$ for @f$ k=1,\dots,K @f$,
     *   - a list of operators @f${\cal L}_{1}, {\cal L}_2, \dots@f$, and
     *   - a list of functions @f$f_1, f_2, \dots@f$,
     *   .
     * then this class helps you computing
     * - gramian matrices, that can look e.g. like the matrices
     *   @f$M_1, M_2, M_3, M_4 \in \mathbb{R}^{K\times K}@f$ with entries
     *   @f{eqnarray*}{
     *    (M_1)_{ij} &=& \langle \varphi_i, \varphi_j \rangle, \\
     *    (M_2)_{ij} &=& \langle {\cal L}_1 \varphi_i, \varphi_j \rangle, \\
     *    (M_3)_{ij} &=& \langle {\cal L}_1 \varphi_i, {\cal L}_1 \varphi_j \rangle, \\
     *    (M_4)_{ij} &=& \langle {\cal L}_2 \varphi_i, {\cal L}_1 \varphi_j \rangle
     *   @f}
     *   for @f$i,j=1,\dots,K@f$,
     * - vectors of scalar products that look e.g. like
     *   @f$b_1, b_2 \in \mathbb{R}^K@f$ with entries
     *   @f{eqnarray*}{
     *    (b_1)_{i} &=& \langle \varphi_i, f_1 \rangle, \\
     *    (b_2)_{i} &=& \langle {\cal L}_1 \varphi_i, f_1 \rangle, \\
     *   @f}
     *   for @f$i=1,\dots,K@f$, and
     * - scalar products like
     *   @f[
     *    \langle f_1, f_2 \rangle \in \mathbb{R}.
     *   @f]
     *
     * In order to use the GramianPipeline the following steps need to be executed:
     *
     * - construct an instance of it and initialize it with the
     *   DiscreteFunctionList:
     *   @code
     *     GramianPipeline pipeline(discFuncList);
     *   @endcode
     * - (optional) register discrete functions or operators and store the returned
     *   handles to these.
     *   @code
     *     FuncHandle hFunc1 = pipeline.registerDiscreteFunction(discFunc1);
     *     OpHandle   hOp1   = pipeline.registerDiscreteOperator(discOp1);
     *     OpHandle   hOp2   = pipeline.registerDiscreteOperator(discOp2);
     *   @endcode
     *   The operator handle for the special identity operator is returned by
     *   @code
     *     OpHandle   hIdOp  = pipeline.getIdentityHandle();
     *   @endcode
     * - add at least one gram matrix, vector or scalar computation
     *   @code
     *     pipeline.addSymmetricGramMatrixComputation(hIdOp, hIdOp, matrix1);
     *     pipeline.addGramMatrixComputation(hOp1, hIdOp, matrix2);
     *     pipeline.addSymmetricGramMatrixComputation(hOp1, hOp1, matrix3);
     *     pipeline.addGramMatrixComputation(hOp2, hOp1, matrix4);
     *     pipeline.addVectorComputation(hIdOp, hFunc1, vector1);
     *     pipeline.addVectorComputation(hOp1, hFunc1, vector2);
     *     pipeline.addScalarComputation(hFunc1, hFunc1, scalar);
     *   @endcode
     * - start the computations of the scalar products
     *   @code
     *     pipeline.run();
     *   @endcode
     *
     *
     * @tparam DiscFuncList   type of the discrete function list
     * @tparam MatrixImp      matrix type where the gramian matrix entries should
     *                        be stored in. Possible types are e.g. MXMatrixWrapper
     *                        or CMatrixWrapper
     */
    template<class DiscFuncList, class MatrixImp>
    class GramianPipeline {
    public:
      //! discrete function list
      typedef DiscFuncList                                           DiscreteFunctionListType;
      //! underlying discrete function
      typedef typename DiscreteFunctionListType
                :: DiscreteFunctionType                              DiscreteFunctionType;
      //! the matrix type for the gramian matrices (c.f. MXMatrixWrapper or CMatrixWrapper)
      typedef MatrixImp                                              MatrixType;
      //! discrete function space
      typedef typename DiscreteFunctionType
                :: DiscreteFunctionSpaceType                         DiscreteFunctionSpaceType;
      //! field type for gramian entries
      typedef typename DiscreteFunctionSpaceType :: RangeFieldType   FieldType;
      //! grid part type
      typedef typename DiscreteFunctionSpaceType :: GridPartType     GridPartType;
      //! entity type for local function evaluations
      typedef typename GridPartType :: GridType
                :: template Codim< 0 > :: Entity                     EntityType;
      //! iterator type for grid walks
      typedef typename DiscreteFunctionSpaceType :: IteratorType     IteratorType;
      //! domain type of DiscreteFunctionType
      typedef typename DiscreteFunctionSpaceType :: DomainType       DomainType;
      //! range type of DiscreteFunctionType
      typedef typename DiscreteFunctionSpaceType :: RangeType        RangeType;
      //! type of func list block
      typedef FuncListBlock< DiscreteFunctionListType >              FuncListBlockType;

    private:
      // only discrete functions which have a localFunction() method are allowed as arguments.
      typedef HasLocalFunction                                       DiscFuncBaseType;
      typedef void EvaluateLocalFuncType(DiscFuncBaseType *&, const EntityType &,
                                         const DomainType &, RangeType &);

      //! stores a local evaluation method from a LocalFuncWrapper (c.f. there for
      //! more details)
      struct EvaluateLocalFuncBase {
        DiscFuncBaseType      *discFunc;
        EvaluateLocalFuncType *evaluateLocalFunc;

        EvaluateLocalFuncBase(DiscFuncBaseType * func,
                              EvaluateLocalFuncType * evLocal)
          : discFunc(func), evaluateLocalFunc(evLocal) {};

        //! evaluates the underlying discrete function locally.
        //! @param en       Entity for which we want to evaluate the local function
        //! @param arg      domain coordinate
        //! @param dest     result of the evaluation
        void evaluateLocal(const EntityType & en, const DomainType & arg, RangeType & dest) {
          (*evaluateLocalFunc)(discFunc, en, arg, dest);
        }
      };

      //! wrapper around a discrete functions that has local functions.
      //! With the static method addInternalFunc, a local function can be
      //! stored in a std::vector of EvaluateLocalFuncBase objects, which
      //! provide a direct evaluation of the localFunction method.
      template<class FuncType>
      struct LocalFuncWrapper {
        static void evaluateLocalWrapper(DiscFuncBaseType *& b,
                                         const EntityType &en,
                                         const DomainType & localX,
                                         RangeType & ret)
        {
          typedef typename FuncType :: LocalFunctionType                 LocalFunctionType;
          const LocalFunctionType lf = static_cast<FuncType &>(*b).localFunction(en);
          lf.evaluate(localX, ret);
        }

        //! adds a wrapper around a discrete function func, and stores the result
        //! at the back of the given vector evLocFuncArray
        //!
        //! @param func             function to be wrapped
        //! @param evLocFuncArray   vector of wrapped functions
        //!
        //! @return new size of vector structure
        static int addInternalFunc(const FuncType* func, std::vector<EvaluateLocalFuncBase> & evLocFuncArray)
        {
          FuncType *func_mutable = const_cast<FuncType *>(func);
          EvaluateLocalFuncBase base(func_mutable, evaluateLocalWrapper);
          evLocFuncArray.push_back(base);
          return evLocFuncArray.size();
        }
      };

      typedef Mapping< FieldType, FieldType, DiscreteFunctionType,
                       DiscreteFunctionType >                        OperatorBaseType;

    public:

      //! handle to an operator that is registered inside a GramianPipeline via
      //! registerDiscreteOperator
      class OpHandle
      {
      public:
        OpHandle(const int opindex, const int dfListSize)
          : opindex_(opindex), dfListSize_(dfListSize) {};

        OpHandle(const OpHandle & op)
          : opindex_(op.opindex_), dfListSize_(op.dfListSize_) {};

        //! returns first index of precomputed operator evaluations
        unsigned int index_begin() const
        {
          return opindex_ * dfListSize_;
        }

        //! returns index behind the of precomputed operator evaluations
        unsigned int index_end() const
        {
          return (opindex_ + 1) * dfListSize_;
        }

        //! returns the internal operator index.
        int op() const
        {
          return opindex_;
        }

      private:
        int opindex_;
        int dfListSize_;
      };

      //! handle to an operator that is registered inside a GramianPipeline via
      //! registerDiscreteFunction
      class FuncHandle
      {
      public:
        FuncHandle(const int funcIndex)
          : funcIndex_(funcIndex) {};

        //! returns the internal function index
        unsigned int index()  const
        {
          return funcIndex_;
        }

      private:
        int funcIndex_;
      };

    private:

      /** helper class that gives access to all information about a single
       * scalar product computation on an entity. Next to a Computation of type
       * ComputationType it stores the matrix indices where the result is added
       * to and the slot indices of the used discrete functions
       */
      template<class ComputationType>
      struct RealComputation
      {
        RealComputation(ComputationType & comp,
                        unsigned int i1 = 0, unsigned int s1 = 0,
                        unsigned int i2 = 0, unsigned int s2 = 0)
          : comp_(&comp), index1(i1), index2(i2), slot1(s1), slot2(s2) {}

        //! returns the computation
        ComputationType & c()
        {
          return *comp_;
        }

      private:
        ComputationType * comp_;
      public:
        // global index in dfList
        unsigned int index1;
        unsigned int index2;
        // local index in block
        unsigned int slot1;
        unsigned int slot2;
      };

      //! class for gram matrix computations where the resulting matrix is symmetric.
      class SymmetricGramMatrixComputation
      {
      public:
        SymmetricGramMatrixComputation(const OpHandle & hOp1, const OpHandle hOp2, MatrixType & matrix)
          : hOp1_(hOp1), hOp2_(hOp2), matrix_(matrix) {};

        //! computes the matrix entries i, j for the entity en
        template<class LocalFunction1, class LocalFunction2, class Entity>
        void compute(const LocalFunction1 & lf1, const LocalFunction2 & lf2, const Entity & en, int i, int j)
        {
          typedef typename Entity :: Geometry                        GeometryType;
          typedef typename LocalFunction1 :: RangeType               RangeType;
          const GeometryType & geo = en.geometry();
          double entityIntegral = 0.0;
          int order_ = 2; //order for the quadrature
          CachingQuadrature<GridPartType, 0> quadrature(en, order_);
          for (unsigned int quadPoint=0; quadPoint != quadrature.nop(); ++quadPoint)
          {
            double integrationElement = geo.integrationElement(quadrature.point(quadPoint));
            double weight = quadrature.weight(quadPoint);

            RangeType ret1;
            RangeType ret2;
            lf1.evaluate( quadrature.point(quadPoint), ret1 );
            lf2.evaluate( quadrature.point(quadPoint), ret2 );
            entityIntegral += integrationElement * weight * ret1 * ret2;
          }
          matrix_.add(i,j,entityIntegral);
          if (i!=j)
            matrix_.add(j,i,entityIntegral);
        }

        const OpHandle& handle1() const {
          return hOp1_;
        }

        const OpHandle& handle2() const {
          return hOp2_;
        }

        //! checks wether the computation depends on the operators with internal
        //! indices i and j
        bool depends(int i, int j)
        {
          return (i == hOp1_.op() && j == hOp2_.op());
        }

      private:
        OpHandle   hOp1_;
        OpHandle   hOp2_;
        MatrixType matrix_;
      };

      //! class for gram matrix computations
      class GramMatrixComputation
      {
      public:
        GramMatrixComputation(const OpHandle & hOp1, const OpHandle hOp2, MatrixType & matrix)
          : hOp1_(hOp1), hOp2_(hOp2), matrix_(matrix) {};

        //! computes the matrix entries i, j for the entity en
        template<class LocalFunction1, class LocalFunction2, class Entity>
        void compute(const LocalFunction1 & lf1, const LocalFunction2 & lf2, const Entity & en, int i, int j)
        {
          typedef typename Entity :: Geometry                        GeometryType;
          const GeometryType & geo = en.geometry();
          typedef typename LocalFunction1 :: RangeType               RangeType;
          double entityIntegral = 0.0;
          int order_ = 2; //order for the quadrature
          CachingQuadrature<GridPartType, 0> quadrature(en, order_);
          for (unsigned int quadPoint=0; quadPoint != quadrature.nop(); ++quadPoint)
            {
              double integrationElement = geo.integrationElement(quadrature.point(quadPoint));
              double weight             = quadrature.weight(quadPoint);

              RangeType ret1;
              RangeType ret2;
              lf1.evaluate( quadrature.point(quadPoint), ret1 );
              lf2.evaluate( quadrature.point(quadPoint), ret2 );
              entityIntegral += integrationElement * weight * ret1 * ret2;
            }
          matrix_.add(i,j,entityIntegral);
        }

        //! checks wether the computation depends on the operators with internal
        //! indices i and j
        bool depends(int i, int j)
        {
          return (i == hOp1_.op() && j == hOp2_.op());
        }

        const OpHandle& handle1() const {
          return hOp1_;
        }

        const OpHandle& handle2() const {
          return hOp2_;
        }

      private:
        OpHandle   hOp1_;
        OpHandle   hOp2_;
        MatrixType matrix_;
      };

      //! class for "gram vector" computations
      class VectorComputation
      {
      public:
        VectorComputation(const OpHandle & hOp1, const EvaluateLocalFuncBase & func,
                          MatrixType & vector, const int fixCol)
          : hOp_(hOp1), func_(func), vector_(vector), fixCol_(fixCol) {};

        //! computes the matrix entries i, fixCol_ for the entity en
        template<class LocalFunction1, class Entity>
        void compute(const LocalFunction1 & lf, const Entity & en, int i)
        {
          typedef typename LocalFunction1 :: RangeType                   RangeType;
          double entityIntegral = 0.0;
          int order_ = 2; //order for the quadrature
          CachingQuadrature<GridPartType, 0> quadrature(en, order_);

          for (unsigned int quadPoint=0; quadPoint != quadrature.nop(); ++quadPoint)
            {
              double integrationElement = en.geometry().integrationElement(quadrature.point(quadPoint));
              double weight = quadrature.weight(quadPoint);

              RangeType ret1;
              RangeType ret2;
              lf.evaluate( quadrature.point(quadPoint), ret1 );
              func_.evaluateLocal( en, coordinate(quadrature.point(quadPoint)), ret2 );
              entityIntegral += integrationElement * weight * ret1 * ret2;
            }
          vector_.add(i,fixCol_,entityIntegral);
        }

        const OpHandle& handle() const {
          return hOp_;
        }

        //! checks wether the computation depends on the operator with internal
        //! index i
        bool depends(int i)
        {
          return (i == hOp_.op());
        }

      private:
        OpHandle              hOp_;
        EvaluateLocalFuncBase func_;
        MatrixType            vector_;
        int                   fixCol_;
      };

      //! class for the computation of a single scalar product between two discrete
      //! functions.
      class ScalarComputation
      {
      public:
        ScalarComputation(const EvaluateLocalFuncBase & func1,
                          const EvaluateLocalFuncBase & func2,
                          MatrixType & scalar,
                          const int fixRow,
                          const int fixCol)
          : func1_(func1), func2_(func2), scalar_(scalar), fixRow_(fixRow), fixCol_(fixCol) {};

        //! computes the matrix entries fixRow_, fixCol_ on entity en
        template<class Entity>
        void compute(const Entity & en)
        {
          double entityIntegral = 0.0;
          int order_ = 2; //order for the quadrature
          CachingQuadrature<GridPartType, 0> quadrature(en, order_);

          for (unsigned int quadPoint=0; quadPoint != quadrature.nop(); ++quadPoint)
            {
              double integrationElement = en.geometry().integrationElement(quadrature.point(quadPoint));
              double weight = quadrature.weight(quadPoint);

              RangeType ret1;
              RangeType ret2;
              func1_.evaluateLocal( en, coordinate(quadrature.point(quadPoint)), ret1 );
              func2_.evaluateLocal( en, coordinate(quadrature.point(quadPoint)), ret2 );
              entityIntegral += integrationElement * weight * ret1 * ret2;
            }
          scalar_.add(fixRow_, fixCol_,entityIntegral);
        }
      private:
        EvaluateLocalFuncBase func1_;
        EvaluateLocalFuncBase func2_;
        MatrixType            scalar_;
        int                   fixRow_;
        int                   fixCol_;
      };

    public:

      /** @brief initializes a GramianPipeline
       *
       * @param dfList              discrete function list for which the gramians
       *                            should be computed.
       * @param parallel_slots      maximum number of discrete functions that are
       *                            held in memory during the computations. If set
       *                            to zero, this is unbounded.  */
      GramianPipeline(DiscreteFunctionListType & dfList, unsigned int parallel_slots=8)
        : dfList_(dfList),
          num_slots_(parallel_slots),
          dfListSize_(dfList_.size()),
          localFuncArray_(),
          localOpArray_(),
          block1_(dfList_, num_slots_),
          block2_(dfList_, num_slots_)
      {
        if(dfListSize_ == 0) {
          DUNE_THROW(InvalidStateException, "GramianPipeline needs to be initialized with non-empty discrete function list");
        }
        if(num_slots_ < 1) {
          DUNE_THROW(InvalidStateException, "Second argument must be at least 1");
        }
      };

      /** @brief destructor */
      ~GramianPipeline()
      {
      }

      /** @brief registers a discrete function to the pipeline and returns a
       * corresponding FuncHandle
       * @param  df       discrete function to register
       * @return discrete function handle */
      template<class DiscreteFunctionImp>
      FuncHandle registerDiscreteFunction(const DiscreteFunctionImp & df)
      {
        // note: the &df statement is legal because a reference inherits the
        // address of the object it points to.
        int size = LocalFuncWrapper<DiscreteFunctionImp> :: addInternalFunc(&df, localFuncArray_);
        FuncHandle hFunc(size-1);

        return hFunc;
      }

      /** @brief registers a discrete operator to the pipeline and returns a
       * corresponding FuncHandle
       * @param  op         discrete operator to register
       * @return discrete operator handle */
      OpHandle registerDiscreteOperator(const OperatorBaseType & op)
      {
        /*    int size = LocalOpWrapper<DiscreteOperatorType> :: addInternalOp(op, localOpArray_);*/
        localOpArray_.push_back(op);
        OpHandle hOp(localOpArray_.size(), dfListSize_);

        assert(hOp.index_end() > 0);

        return hOp;
      }

      /** @brief adds a symmetric gram matrix computation
       *
       * A symmetric gram matrix computation computes a matrix @f$M@f$ with entries
       * @f[
       *   M_{i,j} := \int L_1[\varphi_i] L_2[\varphi_j],
       * @f]
       *
       * where @f$L_1^tL_2@f$ needs to be self-adjoint.
       *
       * @param op1         OpHandle for operator @f$L_1@f$
       * @param op2         OpHandle for operator @f$L_2@f$
       * @param gramian     matrix @f$M@f$ */
      void addSymmetricGramMatrixComputation(const OpHandle & op1, const OpHandle & op2, MatrixType gramian)
      {
        assert(gramian.rows() == dfListSize_ && gramian.cols() == dfListSize_);
        SymmetricGramMatrixComputation gmc(op1, op2, gramian);
        symGramComputations_.push_back(gmc);
      }

      /** @brief returns an operator handle that represents an identy operator
       * @return handle to an identity operator */
      OpHandle getIdentityHandle()
      {
        return OpHandle(0, dfListSize_);
      }

      /** @brief adds a gram matrix computation
       *
       * A gram matrix computation computes a matrix @f$M@f$ with entries
       * @f[
       *   M_{i,j} := \int L_1[\varphi_i] L_2[\varphi_j].
       * @f]
       *
       * @param op1         OpHandle for operator @f$L_1@f$
       * @param op2         OpHandle for operator @f$L_2@f$
       * @param gramian     matrix @f$M@f$ */
      void addGramMatrixComputation(const OpHandle & op1, const OpHandle & op2, MatrixType gramian)
      {
        assert(gramian.rows() == dfListSize_ && gramian.cols() == dfListSize_);
        GramMatrixComputation gmc(op1, op2, gramian);
        gramComputations_.push_back(gmc);
      }

      /** @brief adds a "gram" vector computation
       *
       * A "gram" vector computation computes a matrix @f$M@f$ with entries
       * @f[
       *   M_{i,c} := \int L[\varphi_i] f.
       * @f]
       *
       * @param op          OpHandle for operator @f$L@f$
       * @param func        FuncHandle for discrete function @f$f@f$
       * @param vector      matrix @f$M@f$
       * @param fixColumn   matrix column @f$c@f$ */
      void addVectorComputation(const OpHandle & op, const FuncHandle & func,
                                MatrixType vector, const int fixColumn = 0)
      {
        assert(func.index() < localFuncArray_.size());
        assert(vector.rows() == dfListSize_ && vector.cols() == 1);
        VectorComputation vc(op, localFuncArray_[func.index()], vector, fixColumn);
        vectorComputations_.push_back(vc);

      }

      /** @brief adds a single scalar product computation
       *
       * A scalar computation computes a matrix entry
       * @f[
       *   M_{r,c} := \int f_1 f_2.
       * @f]
       *
       * @param func1       FuncHandle for discrete function @f$f_1@f$
       * @param func2       FuncHandle for discrete function @f$f_2@f$
       * @param scalar      matrix @f$M@f$
       * @param fixRow      matrix row @f$r@f$
       * @param fixColumn   matrix column @f$c@f$ */
      void addScalarComputation(const FuncHandle & func1, const FuncHandle & func2,
                                MatrixType scalar,
                                const int fixRow = 0, const int fixColumn = 0 )
      {
        assert(func1.index() < localFuncArray_.size() && func2.index() < localFuncArray_.size());
        assert(scalar.rows() == 1 && scalar.cols() == 1);
        ScalarComputation sc(localFuncArray_[func1.index()], localFuncArray_[func2.index()], scalar, fixRow, fixColumn);
        scalarComputations_.push_back(sc);
      }

      /** @brief execute the pipeline */
      void run()
      {

        bool include_autarkic_comps = true; // do this only on first gridWalk call
#ifdef GRIDWALK_TEST
        std::vector<std::vector<unsigned int> > gramcombs(dfList_.size(), std::vector<unsigned int>(dfList_.size(), 0));
        std::vector<std::vector<unsigned int> > symgramcombs(dfList_.size(), std::vector<unsigned int>(dfList_.size(), 0));
        std::vector<unsigned int > veccombs(dfList_.size(),  0);
#endif
        std::vector<RealComputation<GramMatrixComputation> >          rGC;
        std::vector<RealComputation<SymmetricGramMatrixComputation> > rsGC;
        std::vector<RealComputation<VectorComputation> >              rVC;

        //! \todo overthink this
        IdOpEvals< FuncListBlockType, DiscreteFunctionType >* idOp
          = new IdOpEvals< FuncListBlockType, DiscreteFunctionType >(2*num_slots_, block1_, block2_);

        opEvals_.push_back(idOp);
        for (unsigned int i=0; i!=localOpArray_.size(); ++i) {
          OpEvals<DiscreteFunctionType>* ptr
            = new OpEvals<DiscreteFunctionType>(2*num_slots_, dfList_.space());
          opEvals_.push_back(ptr);
        }
        block1_.load(0, num_slots_);
        bool block1Moved;
        do
        {    // move first block of functions
          // compute elements belonging only to functions in block 1
          computeBlockContributionsGram(block1_, block1_, rGC, true);
          computeBlockContributionsGramSym(block1_, block1_, rsGC);
          computeBlockContributionsVector(block1_, rVC);
          const unsigned int endOfBlock1_ = block1_.last();

          if (block1_.last()<dfList_.size()-1) {
            block2_.load(endOfBlock1_+1, endOfBlock1_+num_slots_+1);

            do { //move second block of functions
              computeBlockContributionsGram(block1_, block2_, rGC, false);
              computeBlockContributionsGramSym(block1_, block2_, rsGC);
              // before the grid walk compute all operator evaluations that are possible
              // with the current set of functions in block1_ and block2_
              computeOpEvals();
#ifdef GRIDWALK_TEST
              for (unsigned int i = 0; i < rGC.size(); i++) {
                gramcombs[rGC[i].index1][rGC[i].index2]++;
              }
              for (unsigned int i = 0; i < rsGC.size(); i++) {
                symgramcombs[rsGC[i].index1][rsGC[i].index2]++;
                if(rsGC[i].index1!= rsGC[i].index2)
                  symgramcombs[rsGC[i].index2][rsGC[i].index1]++;
              }
              for (unsigned int i = 0; i < rVC.size(); i++) {
                veccombs[rVC[i].index1]++;
              }
#endif
              // do the actual computation here
              gridWalk(rGC, rsGC, rVC, include_autarkic_comps);
            } while (block2_.nextBlock());



            if (block1_.last()+1 != block2_.first())
              block1Moved = block1_.nextBlock();
            else
            {
              block1Moved=false;
              // compute inner elements of second block
              //! \todo maybe these computations can be done above already
              computeBlockContributionsGram(block2_, block2_, rGC, true);
              computeBlockContributionsGramSym(block2_, block2_, rsGC);
              computeBlockContributionsVector(block2_, rVC);
#ifdef GRIDWALK_TEST
              for (unsigned int i = 0; i < rGC.size(); i++) {
                gramcombs[rGC[i].index1][rGC[i].index2]++;
              }
              for (unsigned int i = 0; i < rsGC.size(); i++) {
                symgramcombs[rsGC[i].index1][rsGC[i].index2]++;
                if(rsGC[i].index1!= rsGC[i].index2)
                  symgramcombs[rsGC[i].index2][rsGC[i].index1]++;
              }
              for (unsigned int i = 0; i < rVC.size(); i++) {
                veccombs[rVC[i].index1]++;
              }
#endif
              // before the grid walk compute all operator evaluations that are possible
              // with the current set of functions in block1_ and block2_
              computeOpEvals();
              // do the actual computation here
              gridWalk(rGC, rsGC, rVC, include_autarkic_comps);
            }
          }
          else  {
            block1Moved = false;
            computeOpEvals();
            gridWalk( rGC, rsGC, rVC, include_autarkic_comps );
          }
          //! \todo this has to be closed after the gridwalk itself
        } while (block1Moved);

#ifdef GRIDWALK_TEST
        if(gramComputations_.size() > 0)
        {
          for (unsigned int i = 0; i < gramcombs.size(); i++) {
            for (unsigned int j = 0; j < gramcombs.size(); j++) {
              if (gramcombs[i][j] != gramComputations_.size())
              {
                std::ostringstream oss;
                oss << "Some Gram computations have never been executed:\n";
                for (unsigned int i = 0; i < gramcombs.size(); i++) {
                  for (unsigned int j = 0; j < gramcombs.size(); j++) {
                    oss << " " << gramcombs[i][j] << ", ";
                  }
                  oss << "\n";
                }
                std::cerr << oss.str();
                DUNE_THROW(InvalidStateException, oss.str().c_str());
              }
            }
          }
        }
        if(symGramComputations_.size() > 0)
        {
          for (unsigned int i = 0; i < symgramcombs.size(); i++) {
            for (unsigned int j = 0; j < symgramcombs.size(); j++) {
              if (symgramcombs[i][j] != symGramComputations_.size())
              {
                std::ostringstream oss;
                oss << "Some SymmetricGram computations have never been executed:\n";
                for (unsigned int i = 0; i < symgramcombs.size(); i++) {
                  for (unsigned int j = 0; j < symgramcombs.size(); j++) {
                    oss << " " << symgramcombs[i][j] << ", ";
                  }
                  oss << "\n";
                }
                std::cerr << oss.str();
                DUNE_THROW(InvalidStateException, oss.str().c_str());
              }
            }
          }
        }
        if(vectorComputations_.size() > 0)
        {
          for (unsigned int i = 0; i < veccombs.size(); i++) {
            if (veccombs[i] != vectorComputations_.size())
            {
              std::ostringstream oss;
              oss << "Some Vector Gram computations have never been executed:\n";
              for (unsigned int i = 0; i < veccombs.size(); i++) {
                  oss << " " << veccombs[i] << ", ";
              }
              oss << "\n";
              std::cerr << oss.str();
              DUNE_THROW(InvalidStateException, oss.str().c_str());
            }
          }
        }
#endif
        for (unsigned int i=0; i<opEvals_.size(); ++i) {
          delete opEvals_[i];
        }
      }

    private:

      template< class FuncListBlockType >
      void computeBlockContributionsGram(const FuncListBlockType &block1,
                                         const FuncListBlockType &block2,
                                         std::vector<RealComputation<GramMatrixComputation> > &vec,
                                         const bool &blocksAreEqual ) {
        unsigned int block1size = block1.usedSize();
        unsigned int block2size = block2.usedSize();
        for (unsigned int i=0; i<block1size; ++i) {
          for (unsigned int j=0; j<block2size; ++j) {
            unsigned int cc;
            for( cc = 0 ; cc < gramComputations_.size() ; ++cc )
              {
                // global number of the current function i and j
                const unsigned int ii = block1.mapToGlobal(i);
                const unsigned int jj = block2.mapToGlobal(j);
                RealComputation<GramMatrixComputation> real1(gramComputations_[cc],
                                                             ii, i,
                                                             jj, j);
                vec.push_back(real1);
                if (!blocksAreEqual) {
                RealComputation<GramMatrixComputation> real2(gramComputations_[cc],
                                                             jj, j,
                                                             ii, i);
                vec.push_back(real2);
                }
              }
          } // inner loop (block2 functions)
        } // outer loop (block1 functions)
      }

      template< class FuncListBlockType >
      void computeBlockContributionsGramSym(const FuncListBlockType &block1,
                                            const FuncListBlockType &block2,
                                            std::vector<RealComputation<SymmetricGramMatrixComputation> > &vec ) {
        unsigned int block1size = block1.usedSize();
        unsigned int block2size = block2.usedSize();
        for (unsigned int i=0; i<block1size; ++i) {
          for (unsigned int j=0; j<block2size; ++j) {
            unsigned int cc;
            for( cc = 0 ; cc < symGramComputations_.size() ; ++cc )
              {
                // global number of the current function i and j
                const unsigned int ii = block1.mapToGlobal(i);
                const unsigned int jj = block2.mapToGlobal(j);
                if(ii <= jj)
                  {
                    RealComputation<SymmetricGramMatrixComputation> real1(symGramComputations_[cc],
                                                                        ii, i,
                                                                        jj, j);
                    vec.push_back(real1);
                  }
              }
          } // inner loop (block12functions)
        } // outer loop (block1 functions)
      }

      template< class FuncListBlockType>
      void computeBlockContributionsVector(const FuncListBlockType &block1,
                                           std::vector<RealComputation<VectorComputation> > &vec ) {
        unsigned int block1size = block1.usedSize();
        for (unsigned int i=0; i<block1size; ++i) {
          unsigned int cc;
          for( cc = 0 ; cc < vectorComputations_.size() ; ++cc )
            {
              // global number of the current function i and j
              const unsigned int ii = block1.mapToGlobal(i);
              RealComputation<VectorComputation> real1(vectorComputations_[cc],
                                                       ii, i);
              vec.push_back(real1);
            }
        } // outer loop (block1 functions)
      }

      template < class DiscreteFunctionType >
      class OpEvalsBase {

      public:
        virtual DiscreteFunctionType& getEval(const unsigned int func)=0;

        virtual ~OpEvalsBase() {};

      };

      template< class DiscreteFunctionType >
      class OpEvals : public OpEvalsBase<DiscreteFunctionType> {
        typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      public:
        OpEvals(const unsigned int size, const DiscreteFunctionSpaceType &discFuncSpace)
          :size_(size)
        {
          discFuncSpace_ = &discFuncSpace;
          evals_.resize(size_);
          for (unsigned int i=0; i!=size_; ++i) {
            evals_[i]= new DiscreteFunctionType("opEval", *discFuncSpace_);
          }
        }

        virtual ~OpEvals() {
          for (unsigned int i=0; i!=size_; ++i)
            if (evals_[i]!=NULL)
              delete evals_[i];
        }
      private:
        OpEvals(OpEvals& other) {};

      public:
        void clear()
        {
          map_.clear();
        }

        void addRange(const unsigned int first, const unsigned int last) {
          assert(map_.size()+last-first<size_);
          int j=map_.size();
          for (unsigned int i = first; i<=last; ++i, ++j)
            map_[i]=j;
        }

        /** @brief Get a function
         *
         * @param[in] func global index (index in function list) of desired function.
         *
         * @return returns the discrete function.
         */
        DiscreteFunctionType& getEval(const unsigned int func) {
          assert(map_.find(func)!=map_.end());
          assert(evals_[map_[func]]!=NULL);
          return *evals_[map_[func]];
        }

      private:
        unsigned int                        size_;
        std::vector<DiscreteFunctionType*>  evals_;
        std::map<int, int>                  map_;
        const DiscreteFunctionSpaceType*    discFuncSpace_;
      };



      template < class BlockType, class DiscreteFunctionType >
      class IdOpEvals : public OpEvalsBase<DiscreteFunctionType> {
      public:

        IdOpEvals(const unsigned int &size,
                  BlockType &block1,
                  BlockType &block2) : size_(size),
                                       block1_(block1),
                                       block2_(block2)
        {};

        virtual ~IdOpEvals() {};

        /** @brief Get a function
         *
         * @param[in] func global index (index in function list) of desired function.
         *
         * @return returns the discrete function.
         */
        DiscreteFunctionType& getEval(const unsigned int func) {
          assert((func>=block1_.first() && func <= block1_.last()) ||
                 (func>=block2_.first() && func <= block2_.last()));
          if (func<=block1_.last())
            return block1_.getFunc(func-block1_.first());
          else
            return block2_.getFunc(func-block2_.first());
        }

      private:
        unsigned int size_;
        BlockType &block1_;
        BlockType &block2_;
      };

      /** @brief compute all operator evaluations
       *
       * Evaluate all registered operators on all functions in the
       * current block of functions.
       */
      void computeOpEvals()
      {
        // for each operator
        int end = localOpArray_.size();
        for ( int op = 1; op<=end; ++op)
        {
          // opEvals_[i] is of type OpEvals<...> for i>0!
          // create mapping for operator evaluations
          dynamic_cast<OpEvals<DiscreteFunctionType>*>(opEvals_[op])->clear();
          dynamic_cast<OpEvals<DiscreteFunctionType>*>(opEvals_[op])->addRange(block1_.first(), block1_.last());
          if(block2_.usedSize() > 0)
            dynamic_cast<OpEvals<DiscreteFunctionType>*>(opEvals_[op])->addRange(block2_.first(), block2_.last());
          unsigned int block1size = block1_.usedSize();
          // do the operator evaluations and store them in opEvals_[op] objects
          for ( unsigned int func = 0; func != block1size; ++func ) {
            DiscreteFunctionType &function = opEvals_[op]->getEval(block1_.mapToGlobal(func));
            localOpArray_[op-1](block1_.getFunc(func), function);
          }
          unsigned int block2size = block2_.usedSize();
          for ( unsigned int func = 0; func != block2size; ++func ) {
            DiscreteFunctionType &function =opEvals_[op]->getEval(block2_.mapToGlobal(func));
            localOpArray_[op-1](block2_.getFunc(func), function);
          }
        }
      };

      void gridWalk( std::vector<RealComputation<GramMatrixComputation> >          &rGC,
                     std::vector<RealComputation<SymmetricGramMatrixComputation> > &rsGC,
                     std::vector<RealComputation<VectorComputation> >              &rVC,
                     bool & include_autarkic_comps) {
#ifndef GRIDWALK_TEST
        const IteratorType gridEnd = dfList_.space().gridPart().template end<0>();
        for( IteratorType gridIt = dfList_.space().gridPart().template begin<0>() ;
             gridIt != gridEnd ; ++gridIt )
          {
            EntityType & en = *gridIt;
            typedef typename DiscreteFunctionType :: LocalFunctionType     LocalFunction;
            unsigned int cc;
            for( cc = 0 ; cc < rGC.size() ; ++cc )
              {
                //                 assert(rGC[cc].slot1 < slots_.size() && rGC[cc].slot2 < slots_.size());
                const unsigned int opIndex1 = rGC[cc].c().handle1().op();
                const unsigned int opIndex2 = rGC[cc].c().handle2().op();
                unsigned int row = rGC[cc].index1;
                unsigned int col = rGC[cc].index2;

                const LocalFunction &lf1 = opEvals_[opIndex1]->getEval(row).localFunction(en);
                const LocalFunction &lf2 = opEvals_[opIndex2]->getEval(col).localFunction(en);

                rGC[cc].c().compute(lf1, lf2, en, row, col);
              }

            // compute symmetric operators
            for( cc = 0 ; cc < rsGC.size() ; ++cc )
              {
                const unsigned int opIndex1 = rsGC[cc].c().handle1().op();
                const unsigned int opIndex2 = rsGC[cc].c().handle2().op();
                unsigned int row = rsGC[cc].index1;
                unsigned int col = rsGC[cc].index2;

                const LocalFunction &lf1 = opEvals_[opIndex1]->getEval(row).localFunction(en);
                const LocalFunction &lf2 = opEvals_[opIndex2]->getEval(col).localFunction(en);

                rsGC[cc].c().compute(lf1, lf2, en, row, col);

              }
            for( cc = 0 ; cc < rVC.size() ; ++cc )
              {
                //                assert(rVC[cc].index1 < slots_.size() );
                const unsigned int opIndex = rVC[cc].c().handle().op();
                unsigned int row = rVC[cc].index1;

                const LocalFunction &lf1 = opEvals_[opIndex]->getEval(row).localFunction(en);
                rVC[cc].c().compute(lf1, en, row);
              }
            if( include_autarkic_comps ) {
              for( cc = 0 ; cc < scalarComputations_.size() ; ++cc ) {
                scalarComputations_[cc].compute(en);
              }
            }
          }
#endif
        rGC.clear();
        rsGC.clear();
        rVC.clear();
        include_autarkic_comps = false;
      }

    private:
      //! reference to discrete function list
      DiscreteFunctionListType                    &dfList_;
      //! maxmimum number of slots
      unsigned int                                 num_slots_;
      //! size of discrete function list
      const unsigned int                           dfListSize_;
      //! vector of registered discrete functions
      std::vector<EvaluateLocalFuncBase>           localFuncArray_;
      //! vector of registered discrete operators
      std::vector<OperatorBaseType>                localOpArray_;
      //! set of internal function indices that are in slots
      std::set<unsigned int>                       slot_indices_;
      //! map that maps an internal function index to it's slot number
      std::map<unsigned int, unsigned int>         map_index_to_slot_;
      //! vector of symmetric gram matrix computations
      std::vector<SymmetricGramMatrixComputation>  symGramComputations_;
      //! vector of gram matrix computations
      std::vector<GramMatrixComputation>           gramComputations_;
      //! vector of gram vector computations
      std::vector<VectorComputation>               vectorComputations_;
      //! vector of scalar computations
      std::vector<ScalarComputation>               scalarComputations_;
      //! vector of operator evaluations
      std::vector<OpEvalsBase<DiscreteFunctionType>* >  opEvals_;
      //! Func List Blocks for the underlying func list
      FuncListBlockType      block1_;
      FuncListBlockType      block2_;
    };

  } // end of namespace Dune :: RB

} // end of namespace Dune

#endif  /*DUNE_RB_GRAMIANPIPELINE_HH__*/

/* vim: set sw=2 et: */
