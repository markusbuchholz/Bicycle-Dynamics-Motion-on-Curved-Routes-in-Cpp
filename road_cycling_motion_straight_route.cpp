// Markus Buchholz, 2023
// g++ road_cycling_motion_straight_route.cpp -o t -I/usr/include/python3.8 -lpython3.8

#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <math.h>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

float p = 1.0;  // air density;
float Cd = 1.0; // drag_coeff;
float A = 0.3;  // projected surface
// float u;  // speed - linear velocity
float m = 80.0; // mass
float g = 9.8;  // gravity
float T = 1.0;
float ur = 0.0045;      // rolling resistance
float M = 132.0;        // maximum torque
float D = 10;           // development
float Omega = 2 * M_PI; // meximum pedalling freqency rotation
float b = 0.26;         // braking factor

float dt = 0.01;

//------------------------------------------------------------------
// hydrodynamic drag force

float fd(float u)
{

    return -0.5 * p * Cd * A * u * u;
}

//------------------------------------------------------------------
// gravity force

float fg(float u)
{

    return -m * g;
}

//------------------------------------------------------------------
// friction force

float ff(float u)
{

    return -ur * m * g;
}

//------------------------------------------------------------------
// pedaling force
// https://en.wikipedia.org/wiki/Heaviside_step_function
// equation (9)/(10)

float fp(float u)
{
    return (float)(2 * M_PI * M / D) * (1 - (float)(2 * M_PI * u) / (D * Omega));
}

//------------------------------------------------------------------
// braking force

float fb(float u)
{

    return -b * m * g;
}

//--------------------------------------------------------------------------------
// dot s
float function1(float s, float u)
{
    return s;
}

//--------------------------------------------------------------------------------
// dot u
// pedaling
float function2(float s, float u)
{
    return (fp(u) + fd(u) + ff(u)) / m;
}

//--------------------------------------------------------------------------------
std::tuple<std::vector<float>, std::vector<float>, std::vector<float>> methodRungeKutta1Diff()
{

    std::vector<float> diffEq1;
    std::vector<float> diffEq2;

    std::vector<float> time;

    // init values
    float x1 = 0.0; // s
    float x2 = 5.0; // u
    float t = 0.0;  // init time

    diffEq1.push_back(x1);
    diffEq2.push_back(x2);
    time.push_back(t);

    for (int ii = 0; ii < 10000; ii++)
    {
        t = t + dt;
        float k11 = function1(x1, x2);
        float k12 = function2(x1, x2);

        float k21 = function1(x1 + dt / 2 * k11, x2 + dt / 2 * k12);
        float k22 = function2(x1 + dt / 2 * k11, x2 + dt / 2 * k12);

        float k31 = function1(x1 + dt / 2 * k21, x2 + dt / 2 * k22);
        float k32 = function2(x1 + dt / 2 * k21, x2 + dt / 2 * k22);

        float k41 = function1(x1 + dt * k31, x2 + dt * k32);
        float k42 = function2(x1 + dt * k31, x2 + dt * k32);

        x1 = x1 + dt / 6.0 * (k11 + 2 * k21 + 2 * k31 + k41);
        x2 = x2 + dt / 6.0 * (k12 + 2 * k22 + 2 * k32 + k42);

        diffEq1.push_back(x1);
        diffEq2.push_back(x2);
        time.push_back(t);
    }

    return std::make_tuple(diffEq1, diffEq2, time);
}

//---------------------------------------------------------------------------------------------------------

void plot2D(std::vector<float> xX, std::vector<float> yY)
{
    plt::title("Speed versus time for a rider on a 1000 m flat and straight route ");
    plt::named_plot("speed", xX, yY);
    plt::xlabel("[s]");
    plt::ylabel("[m/s]");
    plt::legend();
    plt::xlabel("[s]");
    plt::ylabel("[m/s]");
    plt::show();
}

//---------------------------------------------------------------------------------------------------------

int main()
{
    auto vdp = methodRungeKutta1Diff();
    auto xX = std::get<2>(vdp); // time
    auto yY = std::get<1>(vdp); // u
    plot2D(xX, yY);
}
