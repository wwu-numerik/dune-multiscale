#ifndef CONSTANTDIFFUSIONMATRIX_HH
#define CONSTANTDIFFUSIONMATRIX_HH

namespace Dune {
namespace Multiscale {

// linear constant diffusion operator that can be filled with given values
// (e.g. with the pre-computed values of a homogenzid matrix)
template< class FieldMatrixImp >
class ConstantDiffusionMatrix
  : public Dune::Fem::Function< Dune::Multiscale::CommonTraits::FunctionSpaceType,
                                ConstantDiffusionMatrix< FieldMatrixImp > >
{
public:
  typedef Dune::Multiscale::CommonTraits::FunctionSpaceType FunctionSpaceType;
  typedef FieldMatrixImp   FieldMatrixType;

private:
  typedef ConstantDiffusionMatrix< FieldMatrixType > ThisType;
  typedef Dune::Fem::Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType::DomainType        DomainType;
  typedef typename FunctionSpaceType::RangeType         RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;

  typedef DomainFieldType TimeType;

public:
  const FieldMatrixType& A_values_;

public:
  inline explicit ConstantDiffusionMatrix(const FieldMatrixType& A_given_values)
    : A_values_(A_given_values)
  {}

  //! (diffusive) flux = A( x , gradient_of_a_function )
  void diffusiveFlux(const DomainType& /*x*/,
                     const JacobianRangeType& gradient,
                     JacobianRangeType& flux) const
  {
      flux[0][0] = A_values_[0][0] * gradient[0][0] + A_values_[0][1] * gradient[0][1];
      flux[0][1] = A_values_[1][0] * gradient[0][0] + A_values_[1][1] * gradient[0][1];
  }

  //! should not be required, since we are in a fully linear setting
  void jacobianDiffusiveFlux(const DomainType&,
                             const JacobianRangeType&,
                             const JacobianRangeType&,
                             JacobianRangeType&) const
  {
      DUNE_THROW(Dune::NotImplemented,"");
  }

};

} //namespace Multiscale {
} //namespace Dune {

#endif // CONSTANTDIFFUSIONMATRIX_HH
