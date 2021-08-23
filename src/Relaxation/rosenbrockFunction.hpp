
#ifndef LSMS_ROSENBROCKFUNCTION_HPP
#define LSMS_ROSENBROCKFUNCTION_HPP


#include <vector>
#include "base_function.hpp"

namespace lsms {

   class RosenbrockFunction : public BaseFunction<double> {

      void evaluate(const std::vector<double> &values,
                    double &result,
                    std::vector<double> &gradient) override;


   };
}


#endif //LSMS_ROSENBROCKFUNCTION_HPP
