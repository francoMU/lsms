
#include <cmath>
#include <iostream>

int main() {


    double Z1 = 23;
    double Z2 = 24;

    double r1[3] = {0.0, 0.0, 0.0};
    double r2[3] = {0.5, 0.5, 0.5};

    double force[3] = {0.0};

    int t = 500;

    double distance{0.0};

    double T[3]{0.0};

    for (int i = -t; i < t+2; i++) {
        for (int j = -t; j < t+1; j++) {
            for (int h = -t; h < t+1; h++) {

                //
                T[0] = i - r2[0];
                T[1] = j - r2[1];
                T[2] = h - r2[2];

                distance = fabs(std::pow(T[0] * T[0] + T[1] * T[1] + T[2] * T[2], 3.0 / 2.0));

                force[0] += Z1 * Z2 * T[0] / distance;
                force[1] += Z1 * Z2 * T[1] / distance;
                force[2] += Z1 * Z2 * T[2] / distance;

                if (i == 0 && j == 0 && h == 0) {
                    continue;
                }

                distance = fabs(std::pow(i * i + j * j + h * h, 3.0 / 2.0));

                force[0] += Z1 * Z1 * i / distance;
                force[1] += Z1 * Z1 * j / distance;
                force[2] += Z1 * Z1 * h / distance;



            }
        }
    }

    std::cout << force[0] << std::endl;
    std::cout << force[1] << std::endl;
    std::cout << force[2] << std::endl;

    //

    double forces_2[3] = {0.0};

    T[0] = 0.2;
    T[1] = 0.5;
    T[2] = 0.5;
    distance = fabs(std::pow(T[0] * T[0] + T[1] * T[1] + T[2] * T[2], 3.0 / 2.0));
    forces_2[0] += Z1 * Z2 * T[0] / distance;
    forces_2[1] += Z1 * Z2 * T[1] / distance;
    forces_2[2] += Z1 * Z2 * T[2] / distance;

    T[0] = 0.2;
    T[1] = -0.5;
    T[2] = 0.5;
    distance = fabs(std::pow(T[0] * T[0] + T[1] * T[1] + T[2] * T[2], 3.0 / 2.0));
    forces_2[0] += Z1 * Z2 * T[0] / distance;
    forces_2[1] += Z1 * Z2 * T[1] / distance;
    forces_2[2] += Z1 * Z2 * T[2] / distance;

    T[0] = 0.2;
    T[1] = -0.5;
    T[2] = -0.5;
    distance = fabs(std::pow(T[0] * T[0] + T[1] * T[1] + T[2] * T[2], 3.0 / 2.0));
    forces_2[0] += Z1 * Z2 * T[0] / distance;
    forces_2[1] += Z1 * Z2 * T[1] / distance;
    forces_2[2] += Z1 * Z2 * T[2] / distance;

    T[0] = 0.2;
    T[1] = 0.5;
    T[2] = -0.5;
    distance = fabs(std::pow(T[0] * T[0] + T[1] * T[1] + T[2] * T[2], 3.0 / 2.0));
    forces_2[0] += Z1 * Z2 * T[0] / distance;
    forces_2[1] += Z1 * Z2 * T[1] / distance;
    forces_2[2] += Z1 * Z2 * T[2] / distance;

    T[0] = -0.8;
    T[1] = 0.5;
    T[2] = 0.5;
    distance = fabs(std::pow(T[0] * T[0] + T[1] * T[1] + T[2] * T[2], 3.0 / 2.0));
    forces_2[0] += Z1 * Z2 * T[0] / distance;
    forces_2[1] += Z1 * Z2 * T[1] / distance;
    forces_2[2] += Z1 * Z2 * T[2] / distance;

    T[0] = -0.8;
    T[1] = -0.5;
    T[2] = 0.5;
    distance = fabs(std::pow(T[0] * T[0] + T[1] * T[1] + T[2] * T[2], 3.0 / 2.0));
    forces_2[0] += Z1 * Z2 * T[0] / distance;
    forces_2[1] += Z1 * Z2 * T[1] / distance;
    forces_2[2] += Z1 * Z2 * T[2] / distance;

    T[0] = -0.8;
    T[1] = -0.5;
    T[2] = -0.5;
    distance = fabs(std::pow(T[0] * T[0] + T[1] * T[1] + T[2] * T[2], 3.0 / 2.0));
    forces_2[0] += Z1 * Z2 * T[0] / distance;
    forces_2[1] += Z1 * Z2 * T[1] / distance;
    forces_2[2] += Z1 * Z2 * T[2] / distance;

    T[0] = -0.8;
    T[1] = 0.5;
    T[2] = -0.5;
    distance = fabs(std::pow(T[0] * T[0] + T[1] * T[1] + T[2] * T[2], 3.0 / 2.0));
    forces_2[0] += Z1 * Z2 * T[0] / distance;
    forces_2[1] += Z1 * Z2 * T[1] / distance;
    forces_2[2] += Z1 * Z2 * T[2] / distance;

    std::cout << forces_2[0] << std::endl;
    std::cout << forces_2[1] << std::endl;
    std::cout << forces_2[2] << std::endl;


    return 0;
}


