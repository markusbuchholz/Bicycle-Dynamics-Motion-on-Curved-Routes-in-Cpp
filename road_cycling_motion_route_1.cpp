// Markus Buchholz, 2023
// g++ road_cycling_motion_route_1.cpp -o t -I/usr/include/python3.8 -lpython3.8

#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <math.h>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

double p = 1.0;  // air density;
double Cd = 1.0; // drag_coeff;
double A = 0.3;  // projected surface
double m = 80.0; // mass
double g = 9.8;  // gravity
double T = 1.0;
double ur = 0.0045;      // rolling resistance
double M = 132.0;        // maximum torque
double D = 10;           // development
double Omega = 2 * M_PI; // meximum pedalling freqency rotation
double b = 0.26;         // braking factor
double Gamma = 0.6;      // threshold pedaling/braking
double fc = 0.4;          // friction cofficient

double dt = 0.01;
int nr_iter = 10000;

// global
std::vector<double> curvature;

//------------------------------------------------------------
struct R
{

    double x;
    double y;
    double z;
};
//------------------------------------------------------------

std::vector<R> cycleRoute()
{
    auto dt_1 = (double)(110.0 / nr_iter);
    std::vector<R> route;
    std::vector<R> velo;

    for (auto ii = -50.0; ii <= 60.0; ii = ii + dt_1)
    {

        R r_i;
        r_i.x = ii;
        r_i.y = std::pow(ii / 10.0, 3);
        r_i.z = 0.0;
        route.push_back(r_i);
    }

    return route;
}
//------------------------------------------------------------

std::vector<R> veloVector()
{
    // dv = dr/dt
    auto dt_1 = (double)(110.0 / nr_iter);
    std::vector<R> rvect;
    std::vector<R> velo;

    for (auto ii = -50.0; ii <= 60.0; ii = ii + dt_1)
    {

        R r_i;
        r_i.x = ii;
        r_i.y = std::pow(ii / 10.0, 3);
        r_i.z = 0.0;
        rvect.push_back(r_i);
    }

    for (int ii = 0; ii < rvect.size() - 1; ii++)
    {
        R v_i;

        v_i.x = (rvect[ii + 1].x - rvect[ii].x) / dt;
        v_i.y = (rvect[ii + 1].y - rvect[ii].y) / dt;
        v_i.z = (rvect[ii + 1].z - rvect[ii].z) / dt;
        velo.push_back(v_i);
    }

    return velo;
}

//------------------------------------------------------------

std::vector<double> computeSpeedVector()
{

    std::vector<R> velo = veloVector();
    std::vector<double> speed;

    for (auto &ii : velo)
    {

        speed.push_back(std::sqrt(std::pow(ii.x, 2) + std::pow(ii.y, 2) + std::pow(ii.z, 2)));
    }

    return speed;
}

//------------------------------------------------------------

std::vector<R> computeTangentVector()
{

    std::vector<R> tangent;

    std::vector<R> veloVec = veloVector();
    std::vector<double> speedVec = computeSpeedVector();

    for (int ii = 0; ii < veloVec.size(); ii++)
    {

        R tang;

        tang.x = veloVec[ii].x / speedVec[ii];
        tang.y = veloVec[ii].y / speedVec[ii];
        tang.z = veloVec[ii].z / speedVec[ii];
        tangent.push_back(tang);
    }

    return tangent;
}

//------------------------------------------------------------

std::vector<double> computeDiffTangentModuleVector()
{

    std::vector<R> tangent = computeTangentVector();
    std::vector<double> tangentDiff;

    for (int ii = 0; ii < tangent.size() - 1; ii++)
    {

        R tang;
        tang.x = tangent[ii + 1].x - tangent[ii].x;
        tang.y = tangent[ii + 1].y - tangent[ii].y;
        tang.z = tangent[ii + 1].z - tangent[ii].z;

        tangentDiff.push_back(std::sqrt(std::pow(tang.x, 2) + std::pow(tang.y, 2) + std::pow(tang.z, 2)));
    }

    return tangentDiff;
}

//------------------------------------------------------------
std::vector<double> computeCurvatureVector()
{
    auto tan_mod = computeDiffTangentModuleVector();
    auto speed = computeSpeedVector();

    std::vector<double> curvVec;

    for (int ii = 0; ii < tan_mod.size(); ii++)
    {

        curvVec.push_back((tan_mod[ii] / speed[ii])*1000);
    }

    return curvVec;
}

//------------------------------------------------------------------
// hydrodynamic drag force

double fd(double u)
{

    return -0.5 * p * Cd * A * u * u;
}

//------------------------------------------------------------------
// gravity force

double fg(double u)
{

    return -m * g;
}

//------------------------------------------------------------------
// friction force

double ff(double u)
{

    return -ur * m * g; 
}

//--------------------------------------------------------------------------------
// swich condition

bool switchCondition(double u, int ii)
{
   
    return curvature[ii] * u * u > Gamma * fc * g ? true : false;
}

//------------------------------------------------------------------
// pedaling force
// https://en.wikipedia.org/wiki/Heaviside_step_function
// equation (9)/(10)
double fp_b(double u, int ii)
{

    if (switchCondition(u, ii) == true)
    {
        std::cout << "here " << "\n";
        return -b * m * g;
    }

    return (double)(2 * M_PI * M / D) * (1 - (double)(2 * M_PI * u) / (D * Omega));
}

//------------------------------------------------------------------
// braking force

double fb(double u)
{

    return -b * m * g;
}

//--------------------------------------------------------------------------------
// dot s
double function1(double s, double u)
{
    return s;
}

//--------------------------------------------------------------------------------
// dot u
// pedaling + braking
double function2(double s, double u, int ii)
{

    return (fp_b(u, ii) + fd(u) + ff(u) /*+ fg(u)*/) / m;
}

//--------------------------------------------------------------------------------
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> methodRungeKutta1Diff()
{

    std::vector<double> diffEq1;
    std::vector<double> diffEq2;

    std::vector<double> time;

    // init values
    double x1 = 0.0; // s
    double x2 = 5.0; // u
    double t = 0.0;  // init time

    diffEq1.push_back(x1);
    diffEq2.push_back(x2);
    time.push_back(t);

    for (int ii = 0; ii < nr_iter - 2; ii++)
    {
        t = t + dt;
        double k11 = function1(x1, x2);
        double k12 = function2(x1, x2, ii);

        double k21 = function1(x1 + dt / 2 * k11, x2 + dt / 2 * k12);
        double k22 = function2(x1 + dt / 2 * k11, x2 + dt / 2 * k12, ii);

        double k31 = function1(x1 + dt / 2 * k21, x2 + dt / 2 * k22);
        double k32 = function2(x1 + dt / 2 * k21, x2 + dt / 2 * k22, ii);

        double k41 = function1(x1 + dt * k31, x2 + dt * k32);
        double k42 = function2(x1 + dt * k31, x2 + dt * k32, ii);

        x1 = x1 + dt / 6.0 * (k11 + 2 * k21 + 2 * k31 + k41);
        x2 = x2 + dt / 6.0 * (k12 + 2 * k22 + 2 * k32 + k42);

        diffEq1.push_back(x1);
        diffEq2.push_back(x2);
        time.push_back(t);
    }

    return std::make_tuple(diffEq1, diffEq2, time);
}

//---------------------------------------------------------------------------------------------------------

void plot2D(std::vector<double> xX, std::vector<double> yY)
{
    plt::title("speed");
    plt::named_plot("speed", xX, yY);
    plt::xlabel("time");
    plt::ylabel("speed");
    plt::legend();
    plt::xlabel("time");
    plt::ylabel("speed");
    plt::show();
}

//---------------------------------------------------------------------------------------------------------

int main()
{

    curvature = computeCurvatureVector();
    int current = 0;
    std::vector<double> time;
    auto value = 0;
    for (auto &ii : curvature)
    {
        time.push_back(value);
        value++;
    }
    plot2D(time, curvature);

   

    auto route = cycleRoute();
    std::vector<double> routeX;
    std::vector<double> routeY;
    std::vector<double> routeZ;

    for (auto &ii : route)
    {
        routeX.push_back(ii.x);
        routeY.push_back(ii.y);
        routeZ.push_back(ii.z);
    }

    plot2D(routeX, routeY);

    auto vdp = methodRungeKutta1Diff();
    auto xX = std::get<2>(vdp); // time
    auto yY = std::get<1>(vdp); // u
   
    plot2D(xX, yY);
}
