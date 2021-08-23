
#include "nl_optimization.hpp"

#include <iostream>
#include <utility>
#include <algorithm>
#include <cmath>


lsms::NLOptimization::NLOptimization(lsms::BaseFunction<double> &base_function,
                                     std::vector<double> starting_value,
                                     int max_iterations,
                                     double tolerance,
                                     double initial_sigma)
      : base_function{base_function},
        max_iterations{max_iterations},
        tolerance{tolerance},
        initial_sigma{initial_sigma},
        starting_value{starting_value} {

   grad_0 = starting_value;
   grad_1 = starting_value;

   x_1 = starting_value;
   x_0 = starting_value;

   r_1 = starting_value;
   r_0 = starting_value;

   x_1 = starting_value;
   x_0 = starting_value;

   d_1 = starting_value;
   d_0 = starting_value;

}

void lsms::NLOptimization::start(const std::vector<double> &x_n_0,
                                 std::vector<double> &x_n_1,
                                 std::vector<double> &grad_n_0,
                                 std::vector<double> &grad_n_1) {

   double results;

   // 1) First step
   base_function.evaluate(x_n_0, results, grad_n_0);

   result_history.emplace_back(results);

   int N = starting_value.size();


   for (int i = 0; i < N; i++) {
      r_0[i] = -grad_n_0[i];
      d_0[i] = -grad_n_0[i];
   };

   // 2) Secant method for finding alpha
   double sigma = initial_sigma;

   for (int i = 0; i < N; i++) {
      x_n_1[i] = x_n_0[i] + sigma * d_0[i];
   };

   base_function.evaluate(x_n_1, results, grad_n_1);

   result_history.emplace_back(results);

   double alpha_upper = 0.0;
   double alpha_lower = 0.0;

   for (int i = 0; i < N; i++) {
      alpha_upper += grad_n_0[i] * d_0[i];
      alpha_lower += grad_n_1[i] * d_0[i] -
                     grad_n_0[i] * d_0[i];
   };

   alpha = -sigma * alpha_upper / alpha_lower;

   // 3) Next step
   for (int i = 0; i < N; i++) {
      x_n_1[i] = x_n_0[i] + alpha * d_0[i];
   };

   // 4) Residual
   base_function.evaluate(x_n_1, results, grad_n_1);

   result_history.emplace_back(results);

   for (int i = 0; i < N; i++) {
      r_1[i] = -grad_n_1[i];
   };

   // 5) Beta
   double beta_up = 0.0;
   double beta_lower = 0.0;

   for (int i = 0; i < N; i++) {
      beta_up += r_1[i] * r_1[i];
      beta_lower += r_0[i] * r_0[i];
   }

   double beta = beta_up / beta_lower;

   // 6) Update congujate
   for (int i = 0; i < N; i++) {
      d_1[i] = r_1[i] + beta * d_0[i];
   }

}


void lsms::NLOptimization::update_step(std::vector<double> &x_n_0,
                                       std::vector<double> &x_n_1) {

   int N = starting_value.size();

   // 7) Shift all values
   for (int i = 0; i < N; i++) {
      x_n_0[i] = x_n_1[i];
      d_0[i] = d_1[i];
      r_0[i] = r_1[i];
   }

}

bool lsms::NLOptimization::calculate_minimum(std::vector<double> &x_result, bool save_steps) {

   if (save_steps) {
      x_history.push_back(starting_value);
   }

   start(x_0, x_1, grad_0, grad_1);

   update_step(x_0, x_1);

   grad_0 = grad_1;

   if (save_steps) {
      x_history.push_back(x_1);
   }

   bool converged = false;

   for (int i = 0; i < max_iterations; i++) {

      iteration(x_0, x_1, grad_0, grad_1);

      if (save_steps) {
         x_history.push_back(x_1);
      }

      if (check_convergence(grad_1)) {
         converged = true;
         x_result = x_1;
         break;
      }

      grad_0 = grad_1;

   }

   return converged;

}

bool lsms::NLOptimization::check_convergence(const std::vector<double> &grad) {

   auto max_val = -1.0;

   for (const double &value : grad) {
      if (std::abs(value) > max_val) {
         max_val = std::abs(value);
      }
   }

   if (max_val > tolerance) {
      return false;
   } else {
      return true;
   }

}


void lsms::NLOptimization::iteration(std::vector<double> &x_n_0,
                                     std::vector<double> &x_n_1,
                                     std::vector<double> &grad_n_0,
                                     std::vector<double> &grad_n_1) {

   double results;
   double sigma;

   int N = starting_value.size();


   // 2) Find alpha

   /*
    * Evaluation of x_n_0 and grad_n_0 was already performed
    */

   sigma = -alpha;

   for (int i = 0; i < N; i++) {
      x_n_1[i] = x_n_0[i] + sigma * d_0[i];
   };

   base_function.evaluate(x_n_1, results, grad_n_1);

   result_history.emplace_back(results);

   double alpha_upper = 0.0;
   double alpha_lower = 0.0;

   for (int i = 0; i < N; i++) {
      alpha_upper += grad_n_0[i] * d_0[i];
      alpha_lower += grad_n_1[i] * d_0[i] -
                     grad_n_0[i] * d_0[i];
   };

   alpha = -sigma * alpha_upper / alpha_lower;

   for (int i = 0; i < N; i++) {
      x_n_1[i] = x_n_0[i] + alpha * d_0[i];
   };

   base_function.evaluate(x_n_1, results, grad_n_1);

   result_history.emplace_back(results);

   for (int i = 0; i < N; i++) {
      r_1[i] = -grad_n_1[i];
   };

   //
   double beta_up = 0.0;
   double beta_lower = 0.0;

   for (int i = 0; i < N; i++) {
      beta_up += r_1[i] * r_1[i];
      beta_lower += r_0[i] * r_0[i];;
   }

   double beta = beta_up / beta_lower;

   for (int i = 0; i < N; i++) {
      d_1[i] = r_1[i] + beta * d_0[i];
   }

   // Shift all values
   for (int i = 0; i < N; i++) {
      x_n_0[i] = x_n_1[i];
      d_0[i] = d_1[i];
      r_0[i] = r_1[i];
   };

}

void lsms::NLOptimization::printRelaxationHistory() const {

   std::printf("Relaxation history\n");

   for (int i = 0; i < result_history.size(); i++) {
      std::printf(" %5d  %20.12f\n", i, result_history[i]);
   }

}
