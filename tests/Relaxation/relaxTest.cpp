
#define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <complex>
#include <cmath>
#include <vector>

#include "nl_optimization.hpp"
#include "rosenbrockFunction.hpp"


namespace TestRelaxation {

    void function(const std::vector<double> &x,
                  double &result,
                  std::vector<double> &gradient) {


        result = (x[0] + 2.0) * (x[0] + 2.0);
        +(x[1] - 2.0) * (x[1] - 2.0);

        gradient[0] = 2.0 * x[0] + 4.0;
        gradient[1] = 2.0 * x[1] - 4.0;

    }


    bool basic_with_generic_function() {

        std::vector<double> grad_0(2, 0.0);
        std::vector<double> grad_1(2, 0.0);

        std::vector<double> x_1(2, 0.0);
        std::vector<double> x_0(2, 0.0);

        std::vector<double> r_1(2, 0.0);
        std::vector<double> r_0(2, 0.0);

        std::vector<double> d_1(2, 0.0);
        std::vector<double> d_0(2, 0.0);

        std::vector<double> beta_0(2, 0.0);
        std::vector<double> beta_1(2, 0.0);

        double result{};

        // 1) First step
        function(x_0, result, grad_0);

        for (int i = 0; i < 2; i++) {
            r_0[i] = -grad_0[i];
            d_0[i] = -grad_0[i];
        };

        // secant method for finding alpha
        double sigma = 0.0001;

        for (int i = 0; i < 2; i++) {
            x_1[i] = x_0[i] + sigma * d_0[i];
        };

        function(x_1, result, grad_1);

        double alpha_upper = 0.0;
        double alpha_lower = 0.0;

        for (int i = 0; i < 2; i++) {
            alpha_upper += grad_0[i] * d_0[i];
            alpha_lower += grad_1[i] * d_0[i] -
                           grad_0[i] * d_0[i];
        };

        double alpha = -sigma * alpha_upper / alpha_lower;

        // 3) Next step
        for (int i = 0; i < 2; i++) {
            x_1[i] = x_0[i] + alpha * d_0[i];
        };

        // 4) Residual
        function(x_1, result, grad_1);
        for (int i = 0; i < 2; i++) {
            r_1[i] = -grad_1[i];
        };

        // 5) Beta
        double beta_up = 0.0;
        double beta_lower = 0.0;

        for (int i = 0; i < 2; i++) {
            beta_up += r_1[i] * r_1[i];
            beta_lower += r_0[i] * r_0[i];
        }

        double beta = beta_up / beta_lower;

        // 6) Update congujate
        for (int i = 0; i < 2; i++) {
            d_1[i] = r_1[i] + beta * d_0[i];
        }

        // 7) Shift all values
        for (int i = 0; i < 2; i++) {
            x_0[i] = x_1[i];
            d_0[i] = d_1[i];
            r_0[i] = r_1[i];
        };


        /*
         * Iteration
         */
        for (int iter = 0; iter < 10000; iter++) {

            // 2) Find alpha
            function(x_0, result, grad_0);


            sigma = -alpha;

            for (int i = 0; i < 2; i++) {
                x_1[i] = x_0[i] + sigma * d_0[i];
            };

            function(x_1, result, grad_1);

            alpha_upper = 0.0;
            alpha_lower = 0.0;

            for (int i = 0; i < 2; i++) {
                alpha_upper += grad_0[i] * d_0[i];
                alpha_lower += grad_1[i] * d_0[i] -
                               grad_0[i] * d_0[i];
            };

            alpha = -sigma * alpha_upper / alpha_lower;

            for (int i = 0; i < 2; i++) {
                x_1[i] = x_0[i] + alpha * d_0[i];
            };

            function(x_1, result, grad_1);

            for (int i = 0; i < 2; i++) {
                r_1[i] = -grad_1[i];
            };

            //
            beta_up = 0.0;
            beta_lower = 0.0;

            for (int i = 0; i < 2; i++) {
                beta_up += r_1[i] * r_1[i];
                beta_lower += r_0[i] * r_0[i];;
            }

            beta = beta_up / beta_lower;

            for (int i = 0; i < 2; i++) {
                d_1[i] = r_1[i] + beta * d_0[i];
            }

            // Shift all values
            for (int i = 0; i < 2; i++) {
                x_0[i] = x_1[i];
                d_0[i] = d_1[i];
                r_0[i] = r_1[i];
            };


        }


        return true;
    }


    bool test_rosenbrock() {

        lsms::RosenbrockFunction function;


        lsms::NLOptimization relaxation(
                function,
                {{0.0, 0.1}});

        // First iteration
        std::vector<double> x_0{{0.0, 0.1}};
        std::vector<double> x_1(2);

        std::vector<double> grad_0(2);
        std::vector<double> grad_1(2);

        relaxation.start(x_0, x_1, grad_0, grad_1);

        relaxation.update_step(x_0, x_1);

        grad_0 = grad_1;

        bool converged = false;

        for (int i = 0; i < 10000; i++) {

            relaxation.iteration(x_0, x_1, grad_0, grad_1);

            std::cout << grad_0[0] << " " << grad_0[1] << std::endl;
            std::cout << x_0[0] << " " << x_0[1] << std::endl;

            if (relaxation.check_convergence(grad_1)) {
                converged = true;
                break;
            }

            grad_0 = grad_1;

        }

        return converged;
    }

    bool test_rosenbrock_2() {

        lsms::RosenbrockFunction function;

        lsms::NLOptimization relaxation(
                function,
                {{0.0, 0.1}},
                10000);

        // First iteration
        std::vector<double> x_result(2);

        auto value = relaxation.calculate_minimum(x_result);

        return value;

    }

}


int main(int argc, char *argv[]) {

    auto test_result =
            TestRelaxation::basic_with_generic_function() &&
            TestRelaxation::test_rosenbrock() &&
            TestRelaxation::test_rosenbrock_2();

    if (test_result) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }

}