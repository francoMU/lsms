
#ifndef MUST_FORCES_PARAMS_HPP
#define MUST_FORCES_PARAMS_HPP


/**
 * @brief class containing all parameters used for the force calculation
 */

namespace lsms {

    constexpr static int DEFAULT_I_Z = 20;
    constexpr static int DEFAULT_I_Y = 20;
    constexpr static int DEFAULT_I_X = 20;
    constexpr static double RADIUS = 300;

    class ForceParameters {

    public:

        /**
         * @param i_x Number of periodic cell in x direction
         * @param i_y Number of periodic cell in y direction
         * @param i_z Number of periodic cell in z direction
         */
        explicit ForceParameters(int i_x = DEFAULT_I_X,
                                 int i_y = DEFAULT_I_Y,
                                 int i_z = DEFAULT_I_Y,
                                 double radius = RADIUS)
                : i_x{i_x}, i_y{i_y}, i_z{i_z}, radius{radius} {};

        int i_x{DEFAULT_I_X};
        int i_y{DEFAULT_I_Y};
        int i_z{DEFAULT_I_Z};
        double radius{RADIUS};

    };

}


#endif //MUST_FORCES_PARAMS_HPP
