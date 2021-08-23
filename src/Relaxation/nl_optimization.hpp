
#ifndef MUST_RELAXATION_HPP
#define MUST_RELAXATION_HPP

#include <vector>

#include "base_function.hpp"

namespace lsms {

   class NLOptimization {

   public:

      NLOptimization(BaseFunction<double> &base_function,
                     std::vector<double> starting_value,
                     int max_iterations = 100,
                     double tolerance = 1e-10,
                     double initial_sigma = 0.00001);

      bool calculate_minimum(std::vector<double> &x_result, bool save_steps = true);

      void iteration(std::vector<double> &x_n_0,
                     std::vector<double> &x_n_1,
                     std::vector<double> &grad_n_0,
                     std::vector<double> &grad_n_1);

      void start(const std::vector<double> &x_n_0,
                 std::vector<double> &x_n_1,
                 std::vector<double> &grad_n_0,
                 std::vector<double> &grad_n_1);

      void update_step(std::vector<double> &x_n_0,
                       std::vector<double> &x_n_1);

      bool check_convergence(const std::vector<double> &grad);

      void printRelaxationHistory() const;

   private:

      BaseFunction<double> &base_function;

      int max_iterations;
      double tolerance;
      double initial_sigma;
      std::vector<double> starting_value;

      std::vector<double> grad_0;
      std::vector<double> grad_1;

      std::vector<double> x_1;
      std::vector<double> x_0;

      std::vector<double> r_1;
      std::vector<double> r_0;

      std::vector<double> d_1;
      std::vector<double> d_0;

      std::vector<std::vector<double>> x_history;

      std::vector<double> result_history;

      double alpha{};


   };

}

#endif //MUST_RELAXATION_HPP
