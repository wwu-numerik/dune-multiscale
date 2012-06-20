#ifndef DUNE_PROBLEM_1_HH
#define DUNE_PROBLEM_1_HH

// Definition of the problem that we want to solve. We define:
// - the right hand side of the problem (realised by the class 'RHSFunction')
// - the exact solution of the problem (realised by the class 'ExactSolution')
// - the boundary condition of the problem (realised by the method 'boundaryTreatment' that sets the boundary-points to
// zero)
// - the matrix/tensor at the left hand side of the problem. For instance the unity matrix if we deal we the laplace
// operator. (the tensor of the problem is realised by the class 'Tensor')

// TODO: Gitterauswahl

// Diese Datei ben�tigt einige der includes und typedefs aus poisson.cc --> nicht ohne weiteres austauschbar
// TODO: besseres Probleminterface �berlegen (mit Klassen o. �.)

// using namespace Dune;

// To explain the implementation of the following classes, we need to know the term 'abstract class'.

// In short: An abstract class is only created to be a 'base class' for a set of other classes (the so called 'derived
// classes'). These derived classes typically share a number of methods that can be subsumed by the abstract class.
// Abstract classes contain at least one virtual method, that means a method of the kind: 'virtual "datatype"
// methodname() = 0;'.
// Examples for virtual methods are: 'virtual void evaluate() = 0;' or 'virtual integer sum() = 0;'.
// Virtual methods can not be used, they will lead to error prompts! Therefor it is impossible to create objects of an
// abstract class. To use such a method nevertheless, the virtual method must be inherited and overwhrighten by an
// equally named method of a derived class.

// First we construct a class 'RHSFunction' that inherits the features/methods of the abstract class 'Function'. In
// order to fill for instance the virtual method 'evaluate' (of //the abstract class 'Function') with life, we need to
// define an equally named method 'evaluate' in 'RHSFunction'.
// Unfortunately, if used very often, virtual methods cause long running times. Therefor, to avoid virtual methods for
// the sake of efficiency, Dune-programmers often make use of the so called Barton-Nackman trick: the derived class
// 'RHSFunction' inherits the (virtual) methods of the class 'Function< FuncSpace, RHSFunction >'. But the class
// template
// 'Function' recieves the template parameter 'RHSFunction', which itself is derived of the class 'Function< FuncSpace,
// RHSFunction >'. The snake seems to bite its tail! But it won't cause any errors. Simplified, the following happens:
// if
// the program wants to execute the method 'evaluate(...)' it is normaly searching in 'Function', there it finds out
// that
// the method is virtual and it has to search in the derived class 'RHSFunction'. Such a process is not very efficient.
// By using the Barton-Nackman trick this process of 'searching' is reversed, the base class 'turns into' the derived
// class and the derived class 'turns into' the base class. Since the program always searches at first in the base class
// (now: 'RHSFunction') it will immediately find a non-virtual method 'evaluate()' that can be executed.
// To realize the Barton-Nackman trick we essentially need the order 'asimp()'. This order must be somewhere in the
// orginal base class ('Function'), before the virtual methods! It is responsible for the reversion of base class and
// derived class.

// !RHSFunction defines the right hand side (RHS) of the governing problem (i.e. it defines 'f').
// The value of the right hand side (i.e. the value of 'f') at 'x' is accessed by the method 'evaluate'. That means 'y
// := f(x)' and 'y' is returned. It is only important that 'RHSFunction' knows the function space ('FuncSpace') that it
// is part from. (f \in FunctionSpace)

class RHSFunction
  : public Function< FuncSpace, RHSFunction >

// we derive a class 'RHSFunction' of the class 'Function< FuncSpace, RHSFunction >'
// Note: 'Function' is only a class-template, 'Function< FuncSpace, RHSFunction >' is a class!

// template parameters of 'Function': - 'FuncSpace' the function space where the Function ('RHSFunction') shall be in,
// and - 'RHSFunction': see Barton-Nackman trick.
{
private:
  typedef Function< FuncSpace, RHSFunction > BaseType;

public:
  RHSFunction(FuncSpace& functionSpace)
    : BaseType(functionSpace)
  {}

  // the following function is defined by 'evaluate':
  // f(x1,x2,x3) = 2 ( (x1) - (x1)^2 ) ( (x2) - (x2)^2 ) +
  // 2 ( (x3) - (x3)^2 ) ( (x2) - (x2)^2 ) +
  // 2 ( (x1) - (x1)^2 ) ( (x3) - (x3)^2 )

  // f(x1,x2,x3) = y

  // y will be returened

  void evaluate(const DomainType& x, RangeType& y) const
  // x = (x1,x2,x3), y = f(x1,x2,x3)
  // Note: in C++, arrays are non-permissible return values. ('array evaluate(...)' causes error prompts)
  // Therefor we need to declare 'evaluate' like above to avoid such problems.
  {
    enum { dim = DomainType::dimension };

    y = 0; // at the beginning we set 'y' to 0. It will be modified during each cycle of the loop below till it is f(x).

    // the following algorithm simply imitates the definition of 'f' from above
    for (int i = 0; i < dim; ++i)
    {
      RangeType tmp = 2.0;
      for (int j = 1; j < dim; ++j)
      {
        const int idx = (i + j) % dim;
        tmp *= x[idx] - SQR(x[idx]);
      }
      y += tmp;
    }
  } // evaluate
};

// ! The exact solution of the problem is defined. This enables us to compare (later on) the numerical solution with
// exact solution in order to calculate the EOC (Estimated Order of Convergence).
class ExactSolution
  : public Function< FuncSpace, ExactSolution >
{
  typedef FuncSpace::DomainFieldType TimeType;
  // essentially: 'DomainFieldType' is the type of an entry of a domain-element.
  // But: it is also used if 'u' (the exact solution) has a time-dependency ('u = u(x,t)').
  // This makes sense since the time-dependency is a one-dimensional element of the 'DomainType' and is therefor also an
  // entry of a domain-element.

public:
  ExactSolution(FuncSpace& f)
    : Function< FuncSpace, ExactSolution >(f) {}

  // The following function (the exakt solution of our problem)
  // ! u(x1,x2,x3) = ( (x1) - (x1)^2 ) * ( (x2) - (x2)^2 ) * ( (x3) - (x3)^2 )
  // is again instantiated by a method 'evaluate'.

  // in case 'u' has NO time-dependency use the following method:
  void evaluate(const DomainType& x, RangeType& y) const
  // x = (x1,x2,x3), y = u(x1,x2,x3)
  {
    y = 1.0;
    for (int i = 0; i < DomainType::dimension; i++)
      y *= ( x[i] - SQR(x[i]) );
  }

  // in case 'u' HAS a time-dependency use the following method:
  void evaluate(const DomainType& x, TimeType time, RangeType& y) const {
    evaluate(x, y);
    // since our exact solution has no time-dependency, we simply refer to the first case (no time-dependency)
  }

  // unfortunately GRAPE requires both cases of the method 'evaluate' to be instantiated
};

// (dummy) implementation of the (for instance diffusion-)tensor/matrix that belongs to the problem:

// since we only deal with the poisson problem, the matrix would be simply the unity-matrix;
// But no matter how we define the following class, the results would be the same. This is due to the fact, that
// 'Tensor' is only a dummy in 'LaplaceFEOp< DiscreteFunctionType, Tensor >' (the sole class that even uses 'Tensor').
// Since the template 'LaplaceFEOp' requieres the class 'Tensor', it is indeed necessary that there existes such a class
// (with its 'evaluate'-methods), but nevertheless it is unimportant HOW it is defined (because it is set to NULL within
// LaplaceFEOp).
class Tensor
  : public Function< FuncSpace, Tensor >
{
public:
  // Constructor for Tensor
  Tensor(FuncSpace& dummytensor)
    : Function< FuncSpace, Tensor >(dummytensor)
  {}

  // instantiate all possible cases of the evaluate-method:

  void evaluate(int i, int j, const DomainType& x, RangeType& y) const {
    y[0] = 12345.67890;
    // NOTE: As already mentioned 'y[0]' could be set to any arbitrary value!
  }

  void evaluate(const DomainType& x, RangeType& y) const {
    y[0] = 12345.67890;
  }

  void evaluate(const DomainType& x, double time, RangeType& y) const {
    y[0] = 12345.67890;
  }
};

// to compute u_h, the dicrete solution of the poisson problem, we need to modify the corresponding linear system of
// equations (LSE).

// Simplified: Let us say that we have to solve

// a_11   a_12   a_13  ...  a_1N         u_1           f_1
// a_21   a_22   a_23  ...  a_2N         u_2           f_2
// .                        .            .             .
// .                        .            .      =      .
// .                        .            .             .
//
// a_N1   a_N2   a_N3  ...  a_NN         u_N           f_N

// with    u_h  =  u_1 \phi_1  +  u_2 \phi_2  +  ... +  u_N \phi_N

// and /phi_i, 1<=i<=N are the base functions of the discrete function space.

// Since we have zero boundary values, u_i must be zero if \phi_i is a base function to a corresponding boundary DOF
// ('Degrees of Freedom').
// (Mathematically: Let  x_i \in \partial \Omega  be a boundary DOF => u_h(x_i) = 0. But  phi_j(x_i) = \delta_(ji) (=>
// phi_i(x_i) = 1) and since the phi_i's are base functions (and therefor linear independent!) we conclude that u_i must
// be zero.)

// Let's fix i (the corresponding \phi_i is a base function to a BOUNDARY DOF).
// By keeping this in mind, we substitute the i-th line of the LSE from above by a_ij^(NEW) := \delta_(ij) and f_i^(NEW)
// := 0.
// Solving the new LSE will lead to solutions u_h with zero boundary values, moreover it will still be a good
// approximation for the exact solution.

// We solve:

// a_11              a_12          a_13      ...      a_1N             u_1             f_1
// a_21              a_22          a_23      ...      a_2N             u_2             f_2
// .                                                                   .              .
// .                                                                   .              .
// .                                                                   .              .
// a_i1^(NEW)     a_i2^(NEW)    a_i3^(NEW)   ...     a_iN^(NEW)        u_i      =      f_i^(NEW)
// .                                                                   .              .
// .                                                                   .              .
// .                                                                   .              .
// a_N1              a_N2          a_N3      ...      a_NN             u_N             f_N

// (Since there is generally more than only one boundary DOF, we typically substitute more than only one line of the
// LSE. Any line substituted is displaced exactly the same way - i.e. a_ij^(NEW) := \delta_(ij) and f_i^(NEW) := 0 )

// NOTE: This possibility of boundary treatment only works if we use OEM-CG-solvers, standard DUNE CG-solvers will most
// likly produce wrong approximations!

// !set the dirichlet points to zero; set dofs to zero
template< class EntityType, class DiscreteFunctionType >
void boundaryTreatment(const EntityType& entity, DiscreteFunctionType& discrete_f) {
  typedef typename DiscreteFunctionType::FunctionSpaceType
  DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::LagrangePointSetType
  LagrangePointSetType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  enum { faceCodim = 1 };

  typedef typename GridPartType::IntersectionIteratorType
  IntersectionIteratorType;

  typedef typename LagrangePointSetType::template Codim< faceCodim >
    ::SubEntityIteratorType
  FaceDofIteratorType;

  // ! get the functionspace of discrete_f??
  const DiscreteFunctionSpaceType& discreteFunctionSpace = discrete_f.space();

  const GridPartType& gridPart = discreteFunctionSpace.gridPart();

// Create an IntersectionIterator
  IntersectionIteratorType it = gridPart.ibegin(entity);
  const IntersectionIteratorType endit = gridPart.iend(entity);
  for ( ; it != endit; ++it)
  {
    if ( !it.boundary() )
      continue;

    LocalFunctionType discrete_f_Local = discrete_f.localFunction(entity);
    const LagrangePointSetType& lagrangePointSet
      = discreteFunctionSpace.lagrangePointSet(entity);

    const int face = it.numberInSelf();
    FaceDofIteratorType faceIt
      = lagrangePointSet.template beginSubEntity< faceCodim >(face);
    const FaceDofIteratorType faceEndIt
      = lagrangePointSet.template endSubEntity< faceCodim >(face);
    for ( ; faceIt != faceEndIt; ++faceIt)
      discrete_f_Local[*faceIt] = 0;
  }
} // boundaryTreatment

#endif // ifndef DUNE_PROBLEM_1_HH
