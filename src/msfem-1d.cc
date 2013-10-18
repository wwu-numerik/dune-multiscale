#include <config.h>
// dune-multiscale
// Copyright Holders: Patrick Henning, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/multiscale/common/main_init.hh>

#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>

#include <dune/common/exceptions.hh>

using namespace std;

// #define EPSILON 0.038776
// #define EPSILON 0.019355
// this conflicted with a local var of the same name in  <dune-istl/dune/istl/solvers.hh>
const double EPSILON = 0.07;

#define NUMBER_OF_STEPS 10

template <typename FunctionType>
double integrate(unsigned int N, double leftB, double rightB, FunctionType& function_to_integrate) {
  if (rightB < leftB) {
    std::stringstream msg;
    msg << "Linke Intervallgrenze " << leftB << " kleiner als rechte Intervallgrenze " << rightB << ". Fehler!"
        << std::endl;
    DUNE_THROW(Dune::InvalidStateException, msg.str());
  }

  if (rightB == leftB) {
    return 0.0;
  }

  double h = (rightB - leftB) / N;

  double integral = 0.0;

  if (function_to_integrate.antiderivative_available) {
    integral =
        function_to_integrate.evaluate_antiderivative(rightB) - function_to_integrate.evaluate_antiderivative(leftB);
  } else {
    for (int i = 0; i < N; i += 1) {
      double x_i = leftB + (i * h);
      integral += h * function_to_integrate.evaluate(x_i);
    }
  }
  return integral;
} // integrate

template <typename FunctionType1, typename FunctionType2>
double error_L2(unsigned int N, double leftB, double rightB, FunctionType1& function1, FunctionType2& function2) {
  if (rightB < leftB) {
    std::stringstream msg;
    msg << "Linke Intervallgrenze " << leftB << " kleiner als rechte Intervallgrenze " << rightB << ". Fehler!"
        << std::endl;
    DUNE_THROW(Dune::InvalidStateException, msg.str());
  }

  if (rightB == leftB) {
    return 0.0;
  }

  double h = (rightB - leftB) / N;

  double error = 0.0;

  for (unsigned int i = 0; i < N; i += 1) {
    double x = leftB + (i * h);
    error += h * (function1.evaluate(x) - function2.evaluate(x)) * (function1.evaluate(x) - function2.evaluate(x));
  }

  return sqrt(error);
} // error_L2

// righ hand side function
class SourceFunction {
public:
  double evaluate(const double /*x*/) {
    // klappt momentan nur fuer f=const
    return 1.0;
  }

  double evaluate_antiderivative(const double x) { return x; }

  double evaluate_antiderivative_antiderivative(const double x) { return 0.5 * x * x; }
};

class A {
public:
  double evaluate(const double y) {
    // double a_in_y = 2.0 + sin( 2.0 * M_PI * y );
    double a_in_y = 1.0 / (2.0 + cos(2.0 * M_PI * y));

    return a_in_y;
  }

  double evaluate_antiderivative_of_inverse(const double y) {
    double val = 2.0 * y;

    val += (1.0 / (2.0 * M_PI)) * sin(2.0 * M_PI * y);
    return val;
  }

  double evaluate_antiderivative_antiderivative_of_inverse(const double y) {
    double val = y * y;

    val -= (1.0 / (4.0 * M_PI * M_PI)) * cos(2.0 * M_PI * y);
    return val;
  }
};

class A_epsilon {
public:
  // evaluate a^{\epsilon}
  double evaluate(const double x) {
    A diffusion_coefficient;
    const double x_new = x / EPSILON;

    double a_eps_in_x = diffusion_coefficient.evaluate(x_new);

    return a_eps_in_x;
  } // evaluate

  // evaluate a^{\epsilon}
  double evaluate_antiderivative_of_inverse(const double x) {
    A diffusion_coefficient;
    const double x_new = x / EPSILON;
    double ret = diffusion_coefficient.evaluate_antiderivative_of_inverse(x_new);

    ret *= EPSILON;
    return ret;
  } // evaluate_antiderivative_of_inverse

  // evaluate a^{\epsilon}
  double evaluate_antiderivative_antiderivative_of_inverse(const double x) {
    A diffusion_coefficient;
    const double x_new = x / EPSILON;
    double ret = diffusion_coefficient.evaluate_antiderivative_antiderivative_of_inverse(x_new);

    ret *= EPSILON * EPSILON;
    return ret;
  } // evaluate_antiderivative_antiderivative_of_inverse
};

class A_MsFEM {
private:
  double* value_vector_;

  // leftB = x_0 and rightB = x_1
  double leftB_, rightB_;

  int N_;

public:
  A_MsFEM(double leftB, double rightB, int N) {
    leftB_ = leftB;
    rightB_ = rightB;
    N_ = N;
    value_vector_ = new double[N];
    for (int i = 0; i < N; i += 1)
      value_vector_[i] = 0.0;
  }

  void set_value(const int i, const double val) {
    if ((i >= N_) || (i < 0)) {
      DUNE_THROW(Dune::InvalidStateException, "Error");
    }

    value_vector_[i] = val;
  } // set_value

  int get_cell_index(const double x) {
    double h = (rightB_ - leftB_) / N_;

    return static_cast<int>((x - leftB_) / h);
  }

  double evaluate(const double x) {
    double h = (rightB_ - leftB_) / N_;
    int i = static_cast<int>((x - leftB_) / h);

    return value_vector_[i];
  }
};

class ZeroFunction {
public:
  double evaluate(const double /*x*/) { return 0.0; }

  double evaluate_antiderivative_of_inverse(const double /*x*/) { return 0.0; }

  double evaluate_antiderivative_antiderivative_of_inverse(const double /*x*/) { return 0.0; }
};

class Exact_Solution {
private:
  // leftB = x_0 and rightB = x_1
  double leftB_, rightB_;

public:
  Exact_Solution(double leftB, double rightB) {
    leftB_ = leftB;
    rightB_ = rightB;
  }

  double evaluate(const double x) {
    SourceFunction f;
    double F_x_1 = f.evaluate_antiderivative(rightB_);
    double F_x_0 = f.evaluate_antiderivative(leftB_);
    double F_x = f.evaluate_antiderivative(x);

    A_epsilon a_eps;

    double ad_a_eps_inverse_x_1 = a_eps.evaluate_antiderivative_of_inverse(rightB_);
    double ad_a_eps_inverse_x_0 = a_eps.evaluate_antiderivative_of_inverse(leftB_);
    double ad_a_eps_inverse_x = a_eps.evaluate_antiderivative_of_inverse(x);

    double ad_ad_a_eps_inverse_x_1 = a_eps.evaluate_antiderivative_antiderivative_of_inverse(rightB_);
    double ad_ad_a_eps_inverse_x_0 = a_eps.evaluate_antiderivative_antiderivative_of_inverse(leftB_);
    double ad_ad_a_eps_inverse_x = a_eps.evaluate_antiderivative_antiderivative_of_inverse(x);

    double a_eps_x_1 = a_eps.evaluate(rightB_);

    double g_x_1 = a_eps_x_1 * (ad_a_eps_inverse_x_1 - ad_a_eps_inverse_x_0);
    double g_x = a_eps_x_1 * (ad_a_eps_inverse_x - ad_a_eps_inverse_x_0);

    double v_x =
        F_x_1 * (ad_a_eps_inverse_x - ad_a_eps_inverse_x_0) - ad_a_eps_inverse_x * F_x + ad_a_eps_inverse_x_0 * F_x_0;

    // + \int_{x_0}^x \int(1/A^eps) f
    v_x += (ad_ad_a_eps_inverse_x - ad_ad_a_eps_inverse_x_0); // in dieser Zeile klappt das nur fuer f=1.
    // (man muss hier eine Stammfunktion von \int(1/A^eps) f kennen.)

    double v_x_1 = F_x_1 * (ad_a_eps_inverse_x_1 - ad_a_eps_inverse_x_0) - ad_a_eps_inverse_x_1 * F_x_1 +
                   ad_a_eps_inverse_x_0 * F_x_0;

    // + \int_{x_0}^x \int(1/A^eps) f
    v_x_1 += (ad_ad_a_eps_inverse_x_1 - ad_ad_a_eps_inverse_x_0); // in dieser Zeile klappt das nur fuer f=1.
    // (man muss hier eine Stammfunktion von \int(1/A^eps) f kennen.)

    return v_x - (g_x * (v_x_1 / g_x_1));
  } // evaluate
};

class Homogenized_Solution {
private:
  // leftB = x_0 and rightB = x_1
  double leftB_, rightB_;

public:
  Homogenized_Solution(double leftB, double rightB) {
    leftB_ = leftB;
    rightB_ = rightB;
  }

  double evaluate(const double x) {
    SourceFunction f;

    double F_x_1 = f.evaluate_antiderivative(rightB_);

    double F_F_x = f.evaluate_antiderivative_antiderivative(x);
    double F_F_x_1 = f.evaluate_antiderivative_antiderivative(rightB_);
    double F_F_x_0 = f.evaluate_antiderivative_antiderivative(leftB_);

    double x_0 = leftB_;
    double x_1 = rightB_;

    A a;

    double a_0 = (1.0 / (a.evaluate_antiderivative_of_inverse(1.0) - a.evaluate_antiderivative_of_inverse(0.0)));

    double v_x = F_x_1 * (x - x_0) - F_F_x + F_F_x_0;

    v_x /= a_0;

    double v_x_1 = F_x_1 * (x_1 - x_0) - F_F_x_1 + F_F_x_0;
    v_x_1 /= a_0;

    double g_x_1 = x_1 - x_0;
    double g_x = x - x_0;

    return v_x - (g_x * (v_x_1 / g_x_1));
  } // evaluate
};

int main(int /*argc*/, char** /*argv[]*/) {
  A a;
  A_epsilon a_eps;

  double left_border = 0.0;
  double right_border = 1.0;

  double a_0 = (1.0 / (a.evaluate_antiderivative_of_inverse(1.0) - a.evaluate_antiderivative_of_inverse(0.0)));
  double a_0_msfem = 0.0;

  double h = (right_border - left_border) / NUMBER_OF_STEPS;

  DSC_LOG_INFO << "EPSILON = " << EPSILON << std::endl;
  DSC_LOG_INFO << "H = " << h << std::endl << std::endl;
  DSC_LOG_INFO << "NUMBER_OF_STEPS = " << NUMBER_OF_STEPS << std::endl << std::endl;

  Exact_Solution u_eps(left_border, right_border);
  Homogenized_Solution u_0(left_border, right_border);
  ZeroFunction zero;

  double error_u_eps_and_u_0 = error_L2(500000, left_border, right_border, u_0, u_eps /*zero*/);

  double u_0_L2_Norm = error_L2(200000, left_border, right_border, u_0, zero);
  // std :: cout << "|| u_0 ||_L2 = " << u_0_L2_Norm << std :: endl;

  double u_eps_L2_Norm = error_L2(500000, left_border, right_border, u_eps, zero);
  // std :: cout << "|| u_eps ||_L2 = " << u_eps_L2_Norm << std :: endl;

  DSC_LOG_INFO << "|| u_eps - u_0 ||_L2 = " << error_u_eps_and_u_0 << std::endl;
  DSC_LOG_INFO << "|| u_eps - u_0 ||_L2 relative = " << error_u_eps_and_u_0 / u_eps_L2_Norm << std::endl << std::endl;

  A_MsFEM a_msfem(left_border, right_border, NUMBER_OF_STEPS);

  for (int i = 0; i < NUMBER_OF_STEPS; i += 1) {
    double x_i = left_border + (i * h);
    double x_i_and_1 = left_border + ((i + 1) * h);

    double msfem_a_i =
        h / (a_eps.evaluate_antiderivative_of_inverse(x_i_and_1) - a_eps.evaluate_antiderivative_of_inverse(x_i));
    // std :: cout << "msfem_a_i = " << msfem_a_i << std :: endl;

    a_msfem.set_value(i, msfem_a_i);

    a_0_msfem += (1.0 / msfem_a_i);
  }

  double error_u_eps_and_u_MsFEM_0 = 0.0;
  double error_u_0_and_u_MsFEM_0 = 0.0;

  double error_u_eps_and_u_MsFEM_eps = 0.0;
  double error_u_0_and_u_MsFEM_eps = 0.0;

  int acc_N = 10000 * NUMBER_OF_STEPS;
  double acc_h = (right_border - left_border) / acc_N;

  double v_x = 0.0;
  double v_x_1 = 0.0;
  double g_x = 0.0;
  double g_x_1 = 0.0;

  SourceFunction f;
  double F_x_1 = f.evaluate_antiderivative(right_border);

  for (int j = 0; j < acc_N; j += 1) {
    double x = left_border + (j * acc_h);

    double F_x = f.evaluate_antiderivative(x);
    double a_msfem_x = a_msfem.evaluate(x);
    v_x_1 += acc_h * (1.0 / a_msfem_x) * (F_x_1 - F_x);

    g_x_1 += acc_h * a_msfem.evaluate(right_border - acc_h) * (1.0 / a_msfem_x);
  }

  double* msfem_solution_vector = new double[NUMBER_OF_STEPS + 1];
  msfem_solution_vector[0] = 0.0;

  double old_value_u_msfem_0_x = 0.0;

  for (int j = 0; j < acc_N; j += 1) {
    double x = left_border + (j * acc_h);

    double F_x = f.evaluate_antiderivative(x);
    double a_msfem_x = a_msfem.evaluate(x);

    v_x += acc_h * (1.0 / a_msfem_x) * (F_x_1 - F_x);
    g_x += acc_h * a_msfem.evaluate(right_border - acc_h) * (1.0 / a_msfem_x);

    double u_msfem_0_x = v_x - (g_x * (v_x_1 / g_x_1));

    msfem_solution_vector[a_msfem.get_cell_index(x) + 1] = u_msfem_0_x;

    double derivative_u_msfem_0_x = (u_msfem_0_x - old_value_u_msfem_0_x) / acc_h;
    old_value_u_msfem_0_x = u_msfem_0_x;

    int i = a_msfem.get_cell_index(x);
    double a_aps_x_i_and_1 = a_eps.evaluate(left_border + ((i + 1) * h));

    double v_i_eps_x = a_aps_x_i_and_1 * (a_eps.evaluate_antiderivative_of_inverse(x) -
                                          a_eps.evaluate_antiderivative_of_inverse(left_border + (i * h)));
    v_i_eps_x -= (x - left_border - (i * h));

    double g_i_eps_x = a_aps_x_i_and_1 * (a_eps.evaluate_antiderivative_of_inverse(x) -
                                          a_eps.evaluate_antiderivative_of_inverse(left_border + (i * h)));

    double g_i_eps_x_i_and_1 =
        a_aps_x_i_and_1 * (a_eps.evaluate_antiderivative_of_inverse(left_border + ((i + 1) * h)) -
                           a_eps.evaluate_antiderivative_of_inverse(left_border + (i * h)));

    double v_i_eps_x_i_and_1 = g_i_eps_x_i_and_1 - h;

    double Q_i_eps = v_i_eps_x - (g_i_eps_x * (v_i_eps_x_i_and_1 / g_i_eps_x_i_and_1));

    double u_msfem_eps_x = u_msfem_0_x + (Q_i_eps * derivative_u_msfem_0_x);

    // std :: cout << "derivative_u_msfem_0_x = " << derivative_u_msfem_0_x << std :: endl;

    // std :: cout << "i = " << a_msfem.get_cell_index(x) << std :: endl;

    // if ( x <= 0.15 ){
    // std :: cout << "u_msfem_0(" << x << ") = " << u_msfem_0_x << std :: endl;}

    error_u_eps_and_u_MsFEM_0 += acc_h * (u_msfem_0_x - u_eps.evaluate(x)) * (u_msfem_0_x - u_eps.evaluate(x));
    error_u_0_and_u_MsFEM_0 += acc_h * (u_msfem_0_x - u_0.evaluate(x)) * (u_msfem_0_x - u_0.evaluate(x));

    error_u_eps_and_u_MsFEM_eps += acc_h * (u_msfem_eps_x - u_eps.evaluate(x)) * (u_msfem_eps_x - u_eps.evaluate(x));
    error_u_0_and_u_MsFEM_eps += acc_h * (u_msfem_eps_x - u_0.evaluate(x)) * (u_msfem_eps_x - u_0.evaluate(x));
  }

  error_u_eps_and_u_MsFEM_0 = sqrt(error_u_eps_and_u_MsFEM_0);
  error_u_0_and_u_MsFEM_0 = sqrt(error_u_0_and_u_MsFEM_0);

  DSC_LOG_INFO << "|| u_eps - u_MsFEM_0 ||_L2 = " << error_u_eps_and_u_MsFEM_0 << std::endl;
  DSC_LOG_INFO << "|| u_eps - u_MsFEM_0 ||_L2 relative = " << error_u_eps_and_u_MsFEM_0 / u_eps_L2_Norm << std::endl
               << std::endl;

  DSC_LOG_INFO << "|| u_0 - u_MsFEM_0 ||_L2 = " << error_u_0_and_u_MsFEM_0 << std::endl;
  DSC_LOG_INFO << "|| u_0 - u_MsFEM_0 ||_L2 relative = " << error_u_0_and_u_MsFEM_0 / u_0_L2_Norm << std::endl
               << std::endl;

  error_u_eps_and_u_MsFEM_eps = sqrt(error_u_eps_and_u_MsFEM_eps);
  error_u_0_and_u_MsFEM_eps = sqrt(error_u_0_and_u_MsFEM_eps);

  DSC_LOG_INFO << "|| u_eps - u_MsFEM_eps ||_L2 = " << error_u_eps_and_u_MsFEM_eps << std::endl;
  DSC_LOG_INFO << "|| u_eps - u_MsFEM_eps ||_L2 relative = " << error_u_eps_and_u_MsFEM_eps / u_eps_L2_Norm << std::endl
               << std::endl;

  DSC_LOG_INFO << "|| u_0 - u_MsFEM_eps ||_L2 = " << error_u_0_and_u_MsFEM_eps << std::endl;
  DSC_LOG_INFO << "|| u_0 - u_MsFEM_eps ||_L2 relative = " << error_u_0_and_u_MsFEM_eps / u_0_L2_Norm << std::endl
               << std::endl;

  a_0_msfem = NUMBER_OF_STEPS * (1.0 / a_0_msfem);

  DSC_LOG_INFO << "A^0 = " << a_0 << std::endl;
  // std :: cout << "A^0_MsFEM = " << a_0_msfem << std :: endl;

  return 0;
} // main
