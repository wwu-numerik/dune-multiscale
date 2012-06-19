
//das alles klappt (mathematisch) nur in 2-D!!! Für Tensoren, die:
// 1. irgendwelche Elliptizitätsbedingungen erfüllen
// 2. A(x,y, \xi) = A(y, \xi) (keine Makroabhängigkeit)
// 3. A( . , \xi) ist 1-periodisch

#ifndef DUNE_HOMOGENIZER_HH
#define DUNE_HOMOGENIZER_HH

// where the quadratures are defined 
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/function/common/function.hh>

//for data output:
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>

//!Wird nicht gebraucht, ist sinnlos:
//#include "../../../problems/EllipticProblems/EllipticCellProblems/cellproblem.hh"


//!NOTE Dieser Homogenizer wäre quasi identisch mit dem HMM Assembler!

namespace Dune 
{


#if 1
// since we need to evaluate A( x, \cdot ) to solve cellproblems (in comparison to A( \cdot, \frac{\cdot}{\epsilon} ) for the global problem), we must transform the orginal tensor to be able to use a standard FEM operator for solving cell problems (otherwise: calling the method evaluate(i,j,x,y) within the matrixassembler would evaluate A^{\epsilon} instead of A(x,\cdot) )
template< class FunctionSpaceImp, class DiffusionImp >
class OnePeriodicDiffusion
		: public Dune::Fem::Function< FunctionSpaceImp, OnePeriodicDiffusion< FunctionSpaceImp, DiffusionImp > >
{
public:
	typedef FunctionSpaceImp FunctionSpaceType;
	typedef DiffusionImp TensorType;

private:
	typedef OnePeriodicDiffusion< FunctionSpaceType, DiffusionImp > ThisType;
	typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
	typedef typename FunctionSpaceType :: DomainType DomainType;
	typedef typename FunctionSpaceType :: RangeType RangeType;
	typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

	typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
	typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

	typedef DomainFieldType TimeType;

	static const int dimDomain = DomainType::dimension;

protected:
	const TensorType &tensor_;

public:
	// Constructor
	inline explicit OnePeriodicDiffusion ( const TensorType &tensor)
		: tensor_( tensor )
	{
	}

	void diffusiveFlux ( const DomainType &y,
						 const JacobianRangeType &direction,
						 JacobianRangeType &flux ) const

	{
		Problem::ModelProblemData model_info;
		const double epsilon = model_info.getEpsilon();

		DomainType new_y;
		for(int i = 0; i < dimDomain; ++i)
			new_y[ i ] = epsilon * y[ i ];

		tensor_.diffusiveFlux( new_y , direction, flux );
	}


	// jacobianDiffusiveFlux = A^{\epsilon}( x , position_gradient ) direction_gradient
	void jacobianDiffusiveFlux ( const DomainType &y,
								 const JacobianRangeType &position_gradient,
								 const JacobianRangeType &direction_gradient,
								 JacobianRangeType &flux ) const
	{

		Problem::ModelProblemData model_info;
		const double epsilon = model_info.getEpsilon();

		DomainType new_y;
		for(int i = 0; i < dimDomain; ++i)
			new_y[ i ] = epsilon * y[ i ];

		tensor_.jacobianDiffusiveFlux( new_y , position_gradient, direction_gradient, flux );

	}


	inline void evaluate ( const int i, const int j,
						   const DomainType &x,
						   const TimeType &time,
						   RangeType &y ) const
	{
		std :: cout << "Do not use this evaluate()-method !!" << std :: endl;
		abort();
		evaluate( i, j, x, y);
	}

	inline void evaluate ( const DomainType &x,
						   RangeType &y ) const
	{
		std :: cout << "Do not use this evaluate()-method !!" << std :: endl;
		abort();
		y = 0;
	}


	inline void evaluate ( const DomainType &x,
						   const TimeType &time,
						   RangeType &y ) const
	{
		std :: cout << "Do not use this evaluate()-method !!" << std :: endl;
		abort();
		y = 0;
	}

};
#endif



} // end namespace 
#endif
