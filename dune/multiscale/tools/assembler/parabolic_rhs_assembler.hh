#ifndef PARABOLIC_RHS_ASSEMBLER_HH
#define PARABOLIC_RHS_ASSEMBLER_HH

namespace Dune {

/**
 * build the right hand side for a discrete parabolic problem that is solved with backward Euler:
     Note that in comparison to the elliptic case:
     1. f needs to be muliplied with the time step size and
     2. the solution u_H^{(k)} of the preceeding time step needs to be added to the right hand side

     Wir haben die rechte Seite = f, wollen darauf aber noch eine discrete Function aufaddieren, die wir discFuncToAdd
     nennen.
     am Ende soll also nicht nur die rechte Seite durch f gebildet werden, sondern durch f+discFuncToAdd
     der Wert wird in discreteFunction gespeichert
     add discreteFunction is an output parameter (kind of return value)
     \todo nothing here is ever instantiated...
 */
template< class DiscreteFunctionImp >
class ParambolicRightHandSideAssembler{
template< int polOrd, class FirstSourceTypeType >
void assembleParabolic(const FirstSourceTypeType& f,
                       const RangeType& time_step_size,
                       const DiscreteFunctionType& u_H_k,
                       DiscreteFunctionType& rhsVector) const {
  // rhsVector is the vector that occurs on the right hand side of the linear system of equations that is to solve
  const DiscreteFunctionSpaceType& discreteFunctionSpace
    = rhsVector.space();

  // set rhs to zero:
  rhsVector.clear();

  const IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    // it* Pointer auf ein Element der Entity
    const GeometryType& geometry = (*it).geometry(); // Referenz auf Geometrie

    LocalFunctionType elementOfRHS = rhsVector.localFunction(*it);   // *it zeigt auf ein bestimmtes Element der
                                                                     // entity
    // hier wird sozusagen ein Pointer von localFunction auf rhs erzeugt. Befinden wir uns auf einer bestimmten
    // entity, so berechnet localFunction alle noetigen Werte und speichert sie (da Pointer) in rhs(aktuelleEntity)
    LocalFunctionType elementOf_u_H_k = u_H_k.localFunction(*it);

    const BaseFunctionSetType baseSet // BaseFunctions leben immer auf Refernzelement!!!
      = discreteFunctionSpace.baseFunctionSet(*it);     // *it Referenz auf eine bestimmtes Element der entity. In der
                                                        // ersten Klasse war das Element fest, deshalb konnte man sich
                                                        // dort Pointer sparen. //loeschen: discreteFunctionSpace
                                                        // statt
                                                        // functionSpace

    const CachingQuadrature< GridPartType, 0 > quadrature(*it, polOrd);   // 0 --> codim 0

    const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade (also die Unbekannten)
    for (int i = 0; i < numDofs; ++i)  // Laufe ueber alle Knoten des entity-elements auf dem wir uns befinden
    {
      // the return values:
      RangeType f_x, phi_x;

      JacobianRangeType gradientPhi;

      const int numQuadraturePoints = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        const double det
          = geometry.integrationElement( quadrature.point(quadraturePoint) );

        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_x':
        f.evaluate(geometry.global( quadrature.point(quadraturePoint) ), f_x);

        // evaluate the current base function at the current quadrature point and save its value in 'phi_x':
        baseSet.evaluate(i, quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;

        // evaluate the gradient of the current base function at the current quadrature point and save its value in
        // 'returnGradient':
        baseSet.jacobian(i, quadrature[quadraturePoint], gradientPhi);
        // Bis jetzt nur der Gradient auf dem Referenzelement!!!!!!! Es muss noch transformiert werden, um den
        // Gradienten auf dem echten Element zu bekommen! Das passiert folgendermassen:

        const FieldMatrix< double, dimension, dimension >& inv
          = geometry.jacobianInverseTransposed( quadrature.point(quadraturePoint) );

        // multiply with transpose of jacobian inverse
        gradientPhi[0] = FMatrixHelp::mult(inv, gradientPhi[0]);

        // value of u_H_k:
        RangeType val_to_add;
        elementOf_u_H_k.evaluate(quadrature, quadraturePoint, val_to_add);

        elementOfRHS[i] += time_step_size * det * quadrature.weight(quadraturePoint) * (f_x * phi_x);

        // add the local value of discFuncToAdd to the right hand side
        elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * val_to_add * phi_x;
      }
    }
  }
}  // end method

template< int polOrd, class FirstSourceTypeType >
void assembleParabolic(const FirstSourceTypeType& f,
                       const TimeType& t,
                       const RangeType& time_step_size,
                       const DiscreteFunctionType& u_H_k,
                       DiscreteFunctionType& rhsVector) const {
  // rhsVector is the vector that occurs on the right hand side of the linear system of equations that is to solve
  // rhs ist der Rueckgabewert der funktion 'assemble'. Hierin wird sozusagen die rechte Seite gespeichert
  const DiscreteFunctionSpaceType& discreteFunctionSpace
    = rhsVector.space();

  // set rhs to zero:
  rhsVector.clear();

  const IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    // it* Pointer auf ein Element der Entity
    const GeometryType& geometry = (*it).geometry(); // Referenz auf Geometrie

    LocalFunctionType elementOfRHS = rhsVector.localFunction(*it);   // *it zeigt auf ein bestimmtes Element der
                                                                     // entity
    // hier wird sozusagen ein Pointer von localFunction auf rhs erzeugt. Befinden wir uns auf einer bestimmten
    // entity, so berechnet localFunction alle noetigen Werte und speichert sie (da Pointer) in rhs(aktuelleEntity)
    const LocalFunctionType elementOf_u_H_k = u_H_k.localFunction(*it);

    const BaseFunctionSetType baseSet // BaseFunctions leben immer auf Refernzelement!!!
      = discreteFunctionSpace.baseFunctionSet(*it);     // *it Referenz auf eine bestimmtes Element der entity. In der
                                                        // ersten Klasse war das Element fest, deshalb konnte man sich
                                                        // dort Pointer sparen. //loeschen: discreteFunctionSpace
                                                        // statt
                                                        // functionSpace

    const CachingQuadrature< GridPartType, 0 > quadrature(*it, polOrd);   // 0 --> codim 0

    const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade (also die Unbekannten)
    for (int i = 0; i < numDofs; ++i)  // Laufe ueber alle Knoten des entity-elements auf dem wir uns befinden
    {
      // the return values:
      RangeType f_x, phi_x;

      JacobianRangeType gradientPhi;

      const int numQuadraturePoints = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        const double det
          = geometry.integrationElement( quadrature.point(quadraturePoint) );

        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_x':
        f.evaluate(geometry.global( quadrature.point(quadraturePoint) ), t, f_x);

        // evaluate the current base function at the current quadrature point and save its value in 'phi_x':
        baseSet.evaluate(i, quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;

        // evaluate the gradient of the current base function at the current quadrature point and save its value in
        // 'returnGradient':
        baseSet.jacobian(i, quadrature[quadraturePoint], gradientPhi);
        // Bis jetzt nur der Gradient auf dem Referenzelement!!!!!!! Es muss noch transformiert werden, um den
        // Gradienten auf dem echten Element zu bekommen! Das passiert folgendermassen:

        const FieldMatrix< double, dimension, dimension >& inv
          = geometry.jacobianInverseTransposed( quadrature.point(quadraturePoint) );

        // multiply with transpose of jacobian inverse
        gradientPhi[0] = FMatrixHelp::mult(inv, gradientPhi[0]);

        // value of u_H_k:
        RangeType val_to_add;
        elementOf_u_H_k.evaluate(quadrature, quadraturePoint, val_to_add);

        elementOfRHS[i] += time_step_size * det * quadrature.weight(quadraturePoint) * (f_x * phi_x);

        // add the local value of discFuncToAdd to the right hand side
        elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * val_to_add * phi_x;
      }
    }
  }
}  // end method

template< int polOrd, class FirstSourceTypeType, class SecondSourceType >
void assembleParabolic(const FirstSourceTypeType& f,
                       const SecondSourceType& G,
                       const RangeType& time_step_size,
                       const DiscreteFunctionType& u_H_k,
                       DiscreteFunctionType& rhsVector) const {
  // rhsVector is the vector that occurs on the right hand side of the linear system of equations that is to solve
  // rhs ist der Rueckgabewert der funktion 'assemble'. Hierin wird sozusagen die rechte Seite gespeichert
  const DiscreteFunctionSpaceType& discreteFunctionSpace
    = rhsVector.space();

  // set rhs to zero:
  rhsVector.clear();

  const IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    // it* Pointer auf ein Element der Entity
    const GeometryType& geometry = (*it).geometry(); // Referenz auf Geometrie

    LocalFunctionType elementOfRHS = rhsVector.localFunction(*it);   // *it zeigt auf ein bestimmtes Element der
                                                                     // entity
    // hier wird sozusagen ein Pointer von localFunction auf rhs erzeugt. Befinden wir uns auf einer bestimmten
    // entity, so berechnet localFunction alle noetigen Werte und speichert sie (da Pointer) in rhs(aktuelleEntity)
    const LocalFunctionType elementOf_u_H_k = u_H_k.localFunction(*it);

    const BaseFunctionSetType baseSet // BaseFunctions leben immer auf Refernzelement!!!
      = discreteFunctionSpace.baseFunctionSet(*it);     // *it Referenz auf eine bestimmtes Element der entity. In der
                                                        // ersten Klasse war das Element fest, deshalb konnte man sich
                                                        // dort Pointer sparen. //loeschen: discreteFunctionSpace
                                                        // statt
                                                        // functionSpace

    const CachingQuadrature< GridPartType, 0 > quadrature(*it, polOrd);   // 0 --> codim 0

    const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade (also die Unbekannten)
    for (int i = 0; i < numDofs; ++i)  // Laufe ueber alle Knoten des entity-elements auf dem wir uns befinden
    {
      // the return values:
      RangeType f_x, phi_x;

      RangeType a[dimension][dimension];

      JacobianRangeType gradientPhi;

      // to save: A \nabla PHI_H
      RangeType G_x[dimension];

      // to save: G \cdot \nabla phi_H
      RangeType val = 0;

      const int numQuadraturePoints = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        const double det
          = geometry.integrationElement( quadrature.point(quadraturePoint) );

        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_x':
        f.evaluate(geometry.global( quadrature.point(quadraturePoint) ), f_x);

        // evaluate the current base function at the current quadrature point and save its value in 'phi_x':
        baseSet.evaluate(i, quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;

        // evaluate the gradient of the current base function at the current quadrature point and save its value in
        // 'returnGradient':
        baseSet.jacobian(i, quadrature[quadraturePoint], gradientPhi);
        // Bis jetzt nur der Gradient auf dem Referenzelement!!!!!!! Es muss noch transformiert werden, um den
        // Gradienten auf dem echten Element zu bekommen! Das passiert folgendermassen:

        const FieldMatrix< double, dimension, dimension >& inv
          = geometry.jacobianInverseTransposed( quadrature.point(quadraturePoint) );

        // multiply with transpose of jacobian inverse
        gradientPhi[0] = FMatrixHelp::mult(inv, gradientPhi[0]);

        // set all entries of G_x to zero to delete old data of a former loop cycle
        for (int k = 0; k < dimension; ++k)
        {
          G_x[k] = 0;
        }

        // the same for val:
        val = 0;

        // evaluate the gradient of Phi_H at the current quadrature point and save its value in 'G_x':
        for (int k = 0; k < dimension; ++k)
        {
          G.evaluate(k, geometry.global( quadrature.point(quadraturePoint) ), G_x[k]);
        }

        for (int k = 0; k < dimension; ++k)
        {
          val += G_x[k] * gradientPhi[0][k];
        }

        // value of u_H_k:
        RangeType val_to_add;
        elementOf_u_H_k.evaluate(quadrature, quadraturePoint, val_to_add);

        elementOfRHS[i] += time_step_size * det * quadrature.weight(quadraturePoint) * (f_x * phi_x);

        // add the local value of discFuncToAdd to the right hand side
        elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * val_to_add * phi_x;

        elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * val;
      }
    }
  }
}  // end method

template< int polOrd, class FirstSourceTypeType, class SecondSourceType >
void assembleParabolic(const FirstSourceTypeType& f,
                       const SecondSourceType& G,
                       const TimeType& t,
                       const RangeType& time_step_size,
                       const DiscreteFunctionType& u_H_k,
                       DiscreteFunctionType& rhsVector) const {
  // rhsVector is the vector that occurs on the right hand side of the linear system of equations that is to solve
  // rhs ist der Rueckgabewert der funktion 'assemble'. Hierin wird sozusagen die rechte Seite gespeichert
  const DiscreteFunctionSpaceType& discreteFunctionSpace
    = rhsVector.space();

  // set rhs to zero:
  rhsVector.clear();

  const IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    // it* Pointer auf ein Element der Entity
    const GeometryType& geometry = (*it).geometry(); // Referenz auf Geometrie

    LocalFunctionType elementOfRHS = rhsVector.localFunction(*it);   // *it zeigt auf ein bestimmtes Element der
                                                                     // entity
    // hier wird sozusagen ein Pointer von localFunction auf rhs erzeugt. Befinden wir uns auf einer bestimmten
    // entity, so berechnet localFunction alle noetigen Werte und speichert sie (da Pointer) in rhs(aktuelleEntity)
    const LocalFunctionType elementOf_u_H_k = u_H_k.localFunction(*it);

    const BaseFunctionSetType baseSet // BaseFunctions leben immer auf Refernzelement!!!
      = discreteFunctionSpace.baseFunctionSet(*it);     // *it Referenz auf eine bestimmtes Element der entity. In der
                                                        // ersten Klasse war das Element fest, deshalb konnte man sich
                                                        // dort Pointer sparen. //loeschen: discreteFunctionSpace
                                                        // statt
                                                        // functionSpace

    const CachingQuadrature< GridPartType, 0 > quadrature(*it, polOrd);   // 0 --> codim 0

    const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade (also die Unbekannten)
    for (int i = 0; i < numDofs; ++i)  // Laufe ueber alle Knoten des entity-elements auf dem wir uns befinden
    {
      // the return values:
      RangeType f_x, phi_x;

      RangeType a[dimension][dimension];

      JacobianRangeType gradientPhi;

      // to save: A \nabla PHI_H
      RangeType G_x[dimension];

      // to save: G \cdot \nabla phi_H
      RangeType val = 0;

      const int numQuadraturePoints = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        const double det
          = geometry.integrationElement( quadrature.point(quadraturePoint) );

        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_x':
        f.evaluate(geometry.global( quadrature.point(quadraturePoint) ), t, f_x);

        // evaluate the current base function at the current quadrature point and save its value in 'phi_x':
        baseSet.evaluate(i, quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;

        // evaluate the gradient of the current base function at the current quadrature point and save its value in
        // 'returnGradient':
        baseSet.jacobian(i, quadrature[quadraturePoint], gradientPhi);
        // Bis jetzt nur der Gradient auf dem Referenzelement!!!!!!! Es muss noch transformiert werden, um den
        // Gradienten auf dem echten Element zu bekommen! Das passiert folgendermassen:

        const FieldMatrix< double, dimension, dimension >& inv
          = geometry.jacobianInverseTransposed( quadrature.point(quadraturePoint) );

        // multiply with transpose of jacobian inverse
        gradientPhi[0] = FMatrixHelp::mult(inv, gradientPhi[0]);

        // set all entries of G_x to zero to delete old data of a former loop cycle
        for (int k = 0; k < dimension; ++k)
        {
          G_x[k] = 0;
        }

        // the same for val:
        val = 0;

        // evaluate the gradient of Phi_H at the current quadrature point and save its value in 'G_x':
        for (int k = 0; k < dimension; ++k)
        {
          G.evaluate(k, geometry.global( quadrature.point(quadraturePoint) ), t, G_x[k]);
        }

        for (int k = 0; k < dimension; ++k)
        {
          val += G_x[k] * gradientPhi[0][k];
        }

        // value of u_H_k:
        RangeType val_to_add;
        elementOf_u_H_k.evaluate(quadrature, quadraturePoint, val_to_add);

        elementOfRHS[i] += time_step_size * det * quadrature.weight(quadraturePoint) * (f_x * phi_x);

        // add the local value of discFuncToAdd to the right hand side
        elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * val_to_add * phi_x;

        elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * val;
      }
    }
  }
}  // end method

template< int polOrd, class FirstSourceType, class InitialValueType >
void assembleParabolic(const FirstSourceType& f,
                       const RangeType& time_step_size,
                       const InitialValueType& v_0,
                       DiscreteFunctionType& rhsVector) const {
  // rhsVector is the vector that occurs on the right hand side of the linear system of equations that is to solve
  // rhs ist der Rueckgabewert der funktion 'assemble'. Hierin wird sozusagen die rechte Seite gespeichert
  const DiscreteFunctionSpaceType& discreteFunctionSpace
    = rhsVector.space();

  // set rhs to zero:
  rhsVector.clear();

  const IteratorType endit = discreteFunctionSpace.end();
  for (IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it)
  {
    // it* Pointer auf ein Element der Entity
    const GeometryType& geometry = (*it).geometry(); // Referenz auf Geometrie

    LocalFunctionType elementOfRHS = rhsVector.localFunction(*it);   // *it zeigt auf ein bestimmtes Element der
                                                                     // entity
    // hier wird sozusagen ein Pointer von localFunction auf rhs erzeugt. Befinden wir uns auf einer bestimmten
    // entity, so berechnet localFunction alle noetigen Werte und speichert sie (da Pointer) in rhs(aktuelleEntity)

    const BaseFunctionSetType baseSet // BaseFunctions leben immer auf Refernzelement!!!
      = discreteFunctionSpace.baseFunctionSet(*it);     // *it Referenz auf eine bestimmtes Element der entity. In der
                                                        // ersten Klasse war das Element fest, deshalb konnte man sich
                                                        // dort Pointer sparen. //loeschen: discreteFunctionSpace
                                                        // statt
                                                        // functionSpace

    const CachingQuadrature< GridPartType, 0 > quadrature(*it, polOrd);   // 0 --> codim 0

    const int numDofs = elementOfRHS.numDofs(); // Dofs = Freiheitsgrade (also die Unbekannten)
    for (int i = 0; i < numDofs; ++i)  // Laufe ueber alle Knoten des entity-elements auf dem wir uns befinden
    {
      // the return values:
      RangeType f_x, v_0_x, phi_x;

      JacobianRangeType gradientPhi;

      const int numQuadraturePoints = quadrature.nop();
      for (int quadraturePoint = 0; quadraturePoint < numQuadraturePoints; ++quadraturePoint)
      {
        const double det
          = geometry.integrationElement( quadrature.point(quadraturePoint) );

        // evaluate the Right Hand Side Function f at the current quadrature point and save its value in 'f_x':
        f.evaluate(geometry.global( quadrature.point(quadraturePoint) ), 0.0, f_x);      // t = 0.0

        // evaluate the Initial Value v_0 at the current quadrature point and save its value in 'v_0_x':
        v_0.evaluate(geometry.global( quadrature.point(quadraturePoint) ), v_0_x);

        // evaluate the current base function at the current quadrature point and save its value in 'phi_x':
        baseSet.evaluate(i, quadrature[quadraturePoint], phi_x);   // i = i'te Basisfunktion;

        // evaluate the gradient of the current base function at the current quadrature point and save its value in
        // 'returnGradient':
        baseSet.jacobian(i, quadrature[quadraturePoint], gradientPhi);
        // Bis jetzt nur der Gradient auf dem Referenzelement!!!!!!! Es muss noch transformiert werden, um den
        // Gradienten auf dem echten Element zu bekommen! Das passiert folgendermassen:

        const FieldMatrix< double, dimension, dimension >& inv
          = geometry.jacobianInverseTransposed( quadrature.point(quadraturePoint) );

        // multiply with transpose of jacobian inverse
        gradientPhi[0] = FMatrixHelp::mult(inv, gradientPhi[0]);

        elementOfRHS[i] += time_step_size * det * quadrature.weight(quadraturePoint) * (f_x * phi_x);

        // add the local initial value to the right hand side
        elementOfRHS[i] += det * quadrature.weight(quadraturePoint) * v_0_x * phi_x;
      }
    }
  }
}  // end method
};
#endif // PARABOLIC_RHS_ASSEMBLER_HH
