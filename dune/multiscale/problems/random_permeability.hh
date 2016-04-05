#ifndef DUNE_MULTISCALE_PROBLEMS_RANDOM_PERMEABILITY_HH
#define DUNE_MULTISCALE_PROBLEMS_RANDOM_PERMEABILITY_HH
#include <exception>
#include <array>
#include <vector>
#include <random>
#include <iostream>
#include <assert.h>
#include <cmath>
#include <mpi.h>
#include <fftw3-mpi.h>

#include <dune/xt/common/exceptions.hh>

/// \brief Class for generating random scalar permeability fields on
/// d-dimensional unit cubes with given correlation function.
///
/// Let \f$ \Omega = [0,1]^d \f$, \f$ x \in \Omega \f$, \f$ k(x) \f$
/// the scalar permeability at \f$ x \f$, \f$ Y(x) = \log(k(x)) \f$,
/// and \f$ R(x,y) = E(Y(x)*Y(y)) \f$ the spacial correlation function.
/// We assume that \f$ R(x,y) = f(x-y) = f(y-x) \f$ is symmetric and depends
/// only on the difference of \f$ x \f$  and \f$ y \f$. Then \f$ Y \f$ is
/// created as linear combination
/// \f$ Y(x) = \sum\limits_{k\in I} \eta_k \, a_k \, \psi_k(x) \f$
/// with index set \f$ I = \left\{ 0,\ldots,2N-1 \right\}^d \f$,
/// \f$ N \f$ : number of segments on \f$ [0,1] \f$ per dimension,
/// \f$ \eta_k \f$ : independent \f$ {\cal N}(0,1) \f$ distributed random
/// numbers,
/// \f$ \psi_k(x) = \exp( i \pi k x ) \f$ :
/// Fourier basis on \f$ [0,2]^d \f$,
/// \f$ a_k \f$ : coefficients chosen to reproduce given correlation function.
///
/// \author jan.mohring@itwm.fraunhofer.de
/// \date   2014
///
/// \tparam dim   space dimension
/// \tparam X     point type //ToDo: with traits of corr
/// \tparam R     coefficient type
/// \tparam COR   class providing correlation via method R operator()(X d)
///               where d=x-y is the difference of two related points
template <int DIM, typename X, typename R, typename COR>
class Permeability {

private:
  typedef std::array<double, 2> complex;
  typedef std::vector<complex> cvec;
  typedef std::array<int, DIM> iarr;

public:
  /// Exception
  class Error : public std::exception {
  public:
    Error(const char* message)
      : _message(message) {}
    const char* what() const throw() { return _message; }

  private:
    const char* _message;
  };

  /// Default constructor
  Permeability() {
    _fft = NULL;
    _ifft = NULL;
  }

  /// Copy constructor
  Permeability(const Permeability&) = delete;

  /// Construct basis from parameters.
  /// \param comm     communicator
  /// \param ldbal    load balancer
  /// \param corr     class providing correlation via method R operator()(X d)
  ///                 where d=x-y is the difference of two related points
  /// \param log2Seg  log2 of number of segments on [0,1] per dimension
  /// \param seed     seed
  /// \param overlap  overlap in domain decomposition (default: 1)
  Permeability(MPI_Comm comm, const COR& corr, int log2Seg, int seed, int overlap = 1, double minimal = 1e-8) {
    _fft = NULL;
    _ifft = NULL;
    init(comm, corr, log2Seg, seed, overlap, minimal);
  };

  /// Construct basis from parameters.
  /// \param comm     communicator
  /// \param ldbal    load balancer
  /// \param corr     class providing correlation via method R operator()(X d)
  ///                 where d=x-y is the difference of two related points
  /// \param log2Seg  log2 of number of segments on [0,1] per dimension
  /// \param seed     seed
  /// \param overlap  overlap in domain decomposition (default: 1)
  /// \param minimal  minimal permeability
  void init(MPI_Comm comm, const COR& corr, int log2Seg, int seed, int overlap = 1, double minimal = 1e-8) {

    // Initialize
    _comm = comm;
    _corr = corr;
    _N = 1 << log2Seg;
    _rand = std::default_random_engine(seed);
    _normal = std::normal_distribution<double>(0, 1);
    _overlap = overlap;
    _minimal = minimal;
    _part = 1;
    if (_fft != NULL) {
      fftw_destroy_plan(_fft);
      fftw_destroy_plan(_ifft);
    };

    ptrdiff_t local_size;
    double h;
    X x;
    int _2N = 2 * _N;
    MPI_Comm_size(_comm, &_nProc);
    MPI_Comm_rank(_comm, &_iProc);
    fftw_mpi_init();
    h = 1.0 / _N;
    setRange();

    // Check input
    assert(log2Seg > 0);

    // Create basis functions in 2D
    if (DIM == 2) {
      local_size = fftw_mpi_local_size_2d(_2N, _2N, _comm, &_n0, &_start);
      assert(_n0 == _2N / _nProc);
      _base = cvec(local_size);
      fftw_complex* base = (fftw_complex*)_base.data();
      _layer = cvec(local_size);
      fftw_complex* layer = (fftw_complex*)_layer.data();
      _fft = fftw_mpi_plan_dft_2d(_2N, _2N, base, base, _comm, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT);
      _ifft =
          fftw_mpi_plan_dft_2d(_2N, _2N, layer, layer, _comm, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
      int p = 0;
      for (int i = 0; i < _n0; ++i) {
        x[0] = (i + _start > _N) ? (_2N - _start - i) * h : (_start + i) * h;
        for (int j = 0; j < _2N; ++j, ++p) {
          x[1] = (j > _N) ? (_2N - j) * h : j * h;
          base[p][0] = _corr(x);
          base[p][1] = 0;
        }
      }
      fftw_execute(_fft);

      double factor = 1.0 / (_2N * _2N);
      for (int i = 0; i < _n0 * _2N; ++i) {
        // assert(base[i]>=0); //TODO: clarify
        base[i][0] = sqrt(factor * fabs(base[i][0]));
        base[i][1] = 0;
      }
    }
    // Create basis functions in 3D
    else if (DIM == 3) {
      local_size = fftw_mpi_local_size_3d(_2N, _2N, _2N, _comm, &_n0, &_start);
      assert(_n0 == _2N / _nProc);
      _base = cvec(local_size);
      fftw_complex* base = (fftw_complex*)_base.data();
      _layer = cvec(local_size);
      fftw_complex* layer = (fftw_complex*)_layer.data();
      _fft =
          fftw_mpi_plan_dft_3d(_2N, _2N, _2N, base, base, _comm, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT);
      _ifft = fftw_mpi_plan_dft_3d(_2N, _2N, _2N, layer, layer, _comm, FFTW_BACKWARD,
                                   FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN);
      int p = 0;
      for (int i = 0; i < _n0; ++i) {
        x[0] = (i + _start > _N) ? (_2N - _start - i) * h : (_start + i) * h;
        for (int j = 0; j < _2N; ++j) {
          x[1] = (j > _N) ? (_2N - j) * h : j * h;
          for (int k = 0; k < _2N; ++k, ++p) {
            x[2] = (k > _N) ? (_2N - k) * h : k * h;
            base[p][0] = _corr(x);
            base[p][1] = 0;
          }
        }
      }
      fftw_execute(_fft);

      double factor = 1.0 / (_2N * _2N * _2N);
      for (int i = 0; i < _n0 * _2N * _2N; ++i) {
        // assert(base[i]>=0); //TODO: clarify
        base[i][0] = sqrt(factor * fabs(base[i][0]));
        base[i][1] = 0;
      }
    } else {
      throw(Error("Only dimensions 2 and 3 are implemented."));
    }
  }

  /// Delete object.
  ~Permeability() {
    fftw_destroy_plan(_fft);
    fftw_destroy_plan(_ifft);
  }

  /// Compute number of processors per dimension
  /// \param   nProc       total number of processors
  /// \param   procPerDim  number of processors per dimension
  /// \returns false, if number of processors is not a power of 2
  static bool partition(int nProc, int* procPerDim) {
    double h = std::log2(nProc);
    int log2P = int(h);
    if (log2P != h)
      return false;
    int nAll = log2P / DIM;
    int nMore = log2P - nAll * DIM;
    for (int i = 0; i < DIM; ++i)
      procPerDim[i] = 1 << (nAll + (i < nMore));
    return true;
  }

private:
  /// Compute coordinates of processor on cartesian processor grid
  /// \param iProc       processor index
  /// \param procPerDim  number of processors per dimension
  /// \param pos         position on processor grid
  static void gridPosition(int iProc, const int* procPerDim, int* pos) {
    int off[DIM];
    off[0] = 1;
    for (int i = 1; i < DIM; ++i) {
      off[i] = off[i - 1] * procPerDim[i - 1];
    }
    for (int i = DIM - 1; i >= 0; --i) {
      pos[i] = iProc / off[i];
      iProc -= pos[i] * off[i];
    }
  }

  /// Compute index range of local part of permeability field
  void setRange() {
    int size = 1;
    int procPerDim[DIM];
    int pos[DIM];
    partition(_nProc, procPerDim);
    gridPosition(_iProc, procPerDim, pos);
    for (int i = 0; i < DIM; ++i) {
      int len = _N / procPerDim[i];
      _iMin[i] = std::max(0, len * pos[i] - _overlap);
      _iMax[i] = std::min(_N, len * (pos[i] + 1) + _overlap);
      _size[i] = _iMax[i] - _iMin[i] + 1;
      size *= _size[i];
    }
    _perm = cvec(size);
  }

  /// Redistribute permeability from layers to blocks.
  void redistribute() {
    MPI_Win win;
    int nCplx = sizeof(complex);
    MPI_Win_create(_layer.data(), nCplx * _layer.size(), nCplx, MPI_INFO_NULL, _comm, &win);
    MPI_Win_fence(MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE, win);
    complex* dest = _perm.data();
    if (DIM == 2) {
      int stepDest = _size[1];
      int stepSource = 2 * _N;
      for (int i = _iMin[0]; i <= _iMax[0]; ++i) {
        int rankSource = i / _n0;
        int source = (i % _n0) * stepSource + _iMin[1];
        MPI_Get(dest, stepDest, MPI_DOUBLE_COMPLEX, rankSource, (MPI_Aint)source, stepDest, MPI_DOUBLE_COMPLEX, win);
        dest += stepDest;
      }
    } else if (DIM == 3) {
      int stepDest = _size[2];
      int stepSource = 4 * _N * _N;
      for (int i = _iMin[0]; i <= _iMax[0]; ++i) {
        int rankSource = i / _n0;
        int source = (i % _n0) * stepSource + 2 * _N * _iMin[1] + _iMin[2];
        for (int j = _iMin[1]; j <= _iMax[1]; j++) {
          MPI_Get(dest, stepDest, MPI_DOUBLE_COMPLEX, rankSource, (MPI_Aint)source, stepDest, MPI_DOUBLE_COMPLEX, win);
          dest += stepDest;
          source += 2 * _N;
        }
      }
    } else {
      throw(Error("Only dimensions 2 and 3 are implemented."));
    }
    MPI_Win_fence(MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOSUCCEED, win);
    MPI_Win_free(&win);
  }

public:
  /// Create random permeability field from basis.
  void create() {
    // toggle between real and imaginary part of permeability field.
    // Recompute only if next real part is requested.
    _part = 1 - _part;
    if (_part == 1) {
      return;
    }
    ptrdiff_t i, n;
    n = (DIM == 2) ? _n0 * 2 * _N : _n0 * 4 * _N * _N;

    // Multiply coefficients with N(0,1)-random numbers and apply IFFT
    for (i = 0; i < n; ++i) {
      _layer[i][0] = _normal(_rand) * _base[i][0];
      _layer[i][1] = _normal(_rand) * _base[i][0];
    }
    fftw_execute(_ifft);

    // k = exp(Y)
    for (i = 0; i < n; ++i) {
      _layer[i][0] = exp(_layer[i][0]);
      _layer[i][1] = exp(_layer[i][1]);
    }

    // Redistribute permeability field
    redistribute();
  }

  /// Evaluate permeability field.
  /// \param x position
  /// \return permeability at x
  R operator()(const X& x) const {
    int cell = 0;
    double t[DIM];
    double k;
    for (int i = 0; i < DIM; ++i) {
      double p = x[i] * _N;
      if (p < _iMin[i] || p > _iMax[i])
        DUNE_THROW(Dune::RangeError, "Coordinate p " << p << " OOB: min " << _iMin[i] << " max " << _iMax[i] << "\n");
      p -= _iMin[i];
      int j = floor(p);
      t[i] = p - j;
      cell = _size[i] * cell + j;
    }
    if (DIM == 2) {
      int c00 = cell;
      int c01 = c00 + 1;
      int c10 = c00 + _size[1];
      int c11 = c10 + 1;
      k = ((1 - t[1]) * _perm[c00][_part] + t[1] * _perm[c01][_part]) * (1 - t[0]) +
          ((1 - t[1]) * _perm[c10][_part] + t[1] * _perm[c11][_part]) * t[0];
    } else if (DIM == 3) {
      int c000 = cell;
      int c001 = c000 + 1;
      int c010 = c000 + _size[2];
      int c011 = c010 + 1;
      int c100 = c000 + _size[1] * _size[2];
      int c101 = c100 + 1;
      int c110 = c100 + _size[2];
      int c111 = c110 + 1;
      k = (((1 - t[2]) * _perm[c000][_part] + t[2] * _perm[c001][_part]) * (1 - t[1]) +
           ((1 - t[2]) * _perm[c010][_part] + t[2] * _perm[c011][_part]) * t[1]) *
              (1 - t[0]) +
          (((1 - t[2]) * _perm[c100][_part] + t[2] * _perm[c101][_part]) * (1 - t[1]) +
           ((1 - t[2]) * _perm[c110][_part] + t[2] * _perm[c111][_part]) * t[1]) *
              t[0];
    } else {
      std::cerr << "Not implemented, yet.\n";
      exit(1);
    }
    return std::max(k, _minimal);
  }

  //--- Members ------------------------------------------------------------
private:
  int _iProc;                               ///< index of processor
  int _nProc;                               ///< total number of processors
  int _overlap;                             ///< overlap of subgrids
  int _N;                                   ///< number of segments on [0,1]
  double _minimal;                          ///< minimal permeability
  ptrdiff_t _n0;                            ///< local num. of segments along 1st dim
  ptrdiff_t _start;                         ///< 1st local index along 1st dim.
  COR _corr;                                ///< correlation as function of x-y
  MPI_Comm _comm;                           ///< MPI-communicator
  fftw_plan _fft;                           ///< plan for fast fourier transform
  fftw_plan _ifft;                          ///< plan for inverse fft
  cvec _base;                               ///< base functions
  cvec _layer;                              ///< local layer of permeability field
  cvec _perm;                               ///< permeability field
  int _part;                                ///< 0: real, 1: imaginary part of _perm
  iarr _iMin;                               ///< minimal global indices of subgrid
  iarr _iMax;                               ///< maximal global indices of subgrid
  iarr _size;                               ///< number of nodes per dim in subgrid
  std::default_random_engine _rand;         ///< random number generator
  std::normal_distribution<double> _normal; ///< normal distribution
};

#endif // DUNE_MULTISCALE_PROBLEMS_RANDOM_
