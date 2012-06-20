#ifndef __CMATRIXWRAPPER_HH__
#define __CMATRIXWRAPPER_HH__

#include <cassert>
#include <cstddef>
#include <iomanip>
#include <cmath>
#include <iostream>

#if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)
 #include <boost/test/unit_test.hpp>
 #include <boost/test/floating_point_comparison.hpp>
#endif // if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)

namespace Dune {
namespace Multiscale {
/** @brief interprets a double array as a matrix with column-wise stored data.
   */
class CMatrixWrapper
{
public:
  /**
     * @brief constructor
     *
     * @param entries_ptr  pointer to (pre-allocated) double array
     * @param rows         number of rows
     * @param cols         number of columns
     */
  CMatrixWrapper(double* entries_ptr,
                 const size_t rows,
                 const size_t cols)
    : entries_ptr_(entries_ptr)
      , rows_(rows)
      , cols_(cols)
  {}

  /**
     * @brief constructor setting data to a given double value
     *
     * @param entries_ptr  pointer to (pre-allocated) double array
     * @param rows         number of rows
     * @param cols         number of columns
     * @param value        the initial value for all entries
     */
  CMatrixWrapper(double* entries_ptr,
                 const size_t rows,
                 const size_t cols,
                 double value)
    : entries_ptr_(entries_ptr)
      , rows_(rows)
      , cols_(cols) {
    clear(value);
  }

  inline void clear(double value = 0.0) {
    size_t i;

    for (i = 0; i < rows_ * cols_; ++i)
    {
      entries_ptr_[i] = value;
    }
  } // clear

  /**
     * @brief add scalar to matrix entry
     *
     * @param i  row position
     * @param j  column position
     * @param val  addent
     */
  void add(const size_t i, const size_t j,
           const double& val) {
    assert(i < rows_ && j < cols_ && entries_ptr_ != NULL);
    assert(i * cols_ + j < rows_ * cols_);
    entries_ptr_[j + i * cols_] += val;
  }

  /**
     * @brief set matrix entry
     *
     * @param i    row position
     * @param j    column position
     * @param val  new value
     */
  void set(const size_t i, const size_t j,
           const double& val) {
    assert(i < rows_ && j < cols_ && entries_ptr_ != NULL);
    assert(i * cols_ + j < rows_ * cols_);
    entries_ptr_[j + i * cols_] = val;
  }

  /**
     * @brief get matrix entry
     *
     * @param i    row position
     * @param j    column position
     * @param val  returned matrix value
     */
  void get(const size_t i, const size_t j,
           double& val) const {
    assert(i < rows_ && j < cols_ && entries_ptr_ != NULL);
    assert(i * cols_ + j < rows_ * cols_);
    val = entries_ptr_[j + i * cols_];
  }

  /**
     * @brief get matrix entry
     *
     * @param i row position
     * @param j column position
     * @return returns a constant reference to the matching entry in the matrix
     */
  const double& operator()(const size_t i, const size_t j) const {
    assert(i < rows_ && j < cols_ && entries_ptr_ != NULL);
    assert(i * cols_ + j < rows_ * cols_);
    return entries_ptr_[j + i * cols_];
  }

  #if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)
  boost::test_tools::predicate_result
  #else
  bool
  #endif // if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)
  compare(const CMatrixWrapper& cmp, double smallEps = 8e-13) {
    if ( ( rows() != cmp.rows() ) || ( cols() != cmp.cols() ) )
    {
      #if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)
      boost::test_tools::predicate_result res(false);
      res.message() << "Different number of rows in matrices: ("
                    << rows() << ", " << cols() << ") != (" << cmp.rows() << ", " << cmp.cols() << ")";
      return res;

      #else // if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)
      return false;

      #endif // if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)
    }
    for (unsigned int i = 0; i < rows_; i++)
    {
      for (unsigned int j = 0; j < cols_; j++)
      {
        double val1;
        double val2;
        get(i, j, val1);
        cmp.get(i, j, val2);
        #if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)
        if (!( boost::test_tools::check_is_close( val1, val2, boost::test_tools::percent_tolerance(0.00001) )
               || ( boost::test_tools::check_is_small(val1, smallEps)
                    && boost::test_tools::check_is_small(val2, smallEps) ) )
            )
        {
          boost::test_tools::predicate_result res(false);
          if ( (cols_ < 500) && (rows_ < 500) )
          {
            std::ostringstream temp1, temp2;
            print(temp1);
            cmp.print(temp2);
            res.message() << "\nthis matrix:\n" << temp1.str() << "\nreference matrix:\n" << temp2.str();
          }
          res.message() << "\nmatrix entries differ at (" << i << ", " << j << "): ["
                        << val1 << " != " << val2 << "]";
          return res;
        }
        #else // if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)
        if (std::fabs(val1 - val2) > 1e-10)
        {
          return false;
        }
        #endif // if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)
      }
    }
    #if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)
    boost::test_tools::predicate_result res(true);
    return res;

    #else // if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)
    return true;

    #endif // if defined (HAVE_BOOST) && defined (BOOST_TEST_MODULE)
  } // compare

  void print(std::ostream& os) const {
    os << "[ ";
    for (unsigned int i = 0; i < rows_; i++)
    {
      os << "  ";
      for (unsigned int j = 0; j < cols_; j++)
      {
        double val;
        get(i, j, val);
        os << std::setw(12) << val;
        if (j != cols_ - 1)
        {
          os << ",";
        }
      }
      if (i != rows_ - 1)
        os << " ;...\n";
      else
        os << " ]\n";
    }
  } // print

  /**
     * @brief number of rows in matrix
     *
     * @return  double
     */
  size_t rows() const {
    return rows_;
  }

  /**
     * @brief number of columns in matrix
     *
     * @return  double
     */
  size_t cols() const {
    return cols_;
  }

  /**@brief Print the matrix to an outstream.
     *
     * @param out an ostream*/
  template< class OutStream >
  void print(OutStream& out) const {
    for (size_t row = 0; row != rows_; ++row)
    {
      out << "|";
      for (size_t col = 0; col != cols_; ++col)
      {
        out << entries_ptr_[row * cols_ + col] << " ";
      }
      out << "|" << std::endl;
    }
    out << std::endl;
  } // print

private:
  double* entries_ptr_;
  size_t rows_;
  size_t cols_;
};
} // end of namespace Dune::RB
} // end of namespace Dune

#endif // __CMATRIXWRAPPER_HH__
