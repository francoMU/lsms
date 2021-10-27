#define OPTIM_ENABLE_ARMA_WRAPPERS

#include <iostream>

#include "optim.hpp"


class Data {
public:
  double value{5};
  double *ptr_value;
};

class AckleyFunction {

public:

  static double eval(const arma::vec &vals_inp, arma::vec *grad_out, void *opt_data) {
    const double x = vals_inp(0);
    const double y = vals_inp(1);
    const double pi = arma::datum::pi;

    double obj_val = -20 * std::exp(-0.2 * std::sqrt(0.5 * (x * x + y * y))) -
                     std::exp(0.5 * (std::cos(2 * pi * x) + std::cos(2 * pi * y))) + 22.718282L;


    Data *opt_data_ptr = (Data *) opt_data;

    //std::cout << opt_data_ptr->value << std::endl;
    //std::cout << *(opt_data_ptr->ptr_value) << std::endl;

    return obj_val;
  }

};


int main() {
  arma::vec x = arma::ones(2, 1) + 1.0; // (2,2)


  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

  double value = 10;

  Data data{};
  data.ptr_value = &value;

  bool success = optim::pso(x, AckleyFunction::eval, &data);

  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;

  if (success) {
    std::cout << "de: Ackley test completed successfully.\n"
              << "elapsed time: " << elapsed_seconds.count() << "s\n";
  } else {
    std::cout << "de: Ackley test completed unsuccessfully." << std::endl;
  }

  arma::cout << "\nde: solution to Ackley test:\n" << x << arma::endl;

  return 0;
}