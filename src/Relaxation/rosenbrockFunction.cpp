/**
 *
 * RadialDFTLib
 *
 * rosenbrockFunction.cpp
 *
 * Created by Franco Moitzi on 8/11/21.
 *
 * Copyright (c) 2021 University of Leoben. All rights reserved.
 *
 */
#include "rosenbrockFunction.hpp"


static void eval_value(const std::vector<double> &values,
                       double &result) {

   result = 100 * (values[1] - values[0] * values[0]) + (1 - values[0]) * (1 - values[0]);

}

static void eval_gradient(const std::vector<double> &values,
                          std::vector<double> &result) {

   result[0] = 400 * values[0] * values[0] * values[0]
               - 400 * values[1] * values[0]
               + 2 * values[0]
               - 2;

   result[1] = 200 * (values[1] - values[0] * values[0]);

}

void
lsms::RosenbrockFunction::evaluate(const std::vector<double> &values,
                                   double &result,
                                   std::vector<double> &gradient) {

   eval_gradient(values, gradient);

   eval_value(values, result);
}
