// population dynamics simulation implementations

// used headers and libraries
#include <iostream>
#include "NumericalIntegration.hh"
#include "Vector2D.hh"
#include <vector>

// ---------------------------------------------------------------------------------------------------------

// initial time
double const t0 = 0.;
// integration time boundary
double const t1 = 30.;
// step size
double const h = 1e-2;
// set value for r (will be reset)
double r = 0;

// ---------------------------------------------------------------------------------------------------------

// differences of birth and death for coupled logistic equations
double const r1 = 0.5;
double const r2 = 0.5;
// coupling coefficients for coupled logistic equations (will be reset)
double alpha = 1.;
double beta = 1.;
// saturation values for coupled logistic equations (will be reset)
double k1 = 1.;
double k2 = 1.;

// ---------------------------------------------------------------------------------------------------------

// parameters for the Lotka-Volterra model
double const a = 20.;
double const b = 1.;
double const c = 30.;
double const d = 1.;

// ---------------------------------------------------------------------------------------------------------

// RHS for logistic equation
auto LogisticEq = [&](double, double x) {
    return r * x * (1 - x);
};

// ---------------------------------------------------------------------------------------------------------

// RHS for coupled logistic equations
auto LogisticCoupledEq = [&](double, vector2<double> state) -> vector2<double> {
    return {
        r1 * state.x * (1 - (state.x + alpha * state.y) / k1),
        r2 * state.y * (1 - (state.y + beta * state.x) / k2)};
};

// ---------------------------------------------------------------------------------------------------------

// RHS for Lotka-Volterra model
auto LotkaVolterraEq = [&](double, vector2<double> state) -> vector2<double> {
    return {
        a * state.x - b * state.x * state.y,
        c * state.x * state.y - d * state.y};
};

// ---------------------------------------------------------------------------------------------------------

// callback function for 1D: write to file
auto toFile1D = [&](double t, double x, std::string fileName) {
    std::ofstream data;
    data.open(fileName, std::ofstream::app);
    data << t << " " << x << std::endl;
    data.close();
};

// ---------------------------------------------------------------------------------------------------------

// callback function for 2D: write to file
auto toFile2D = [&](double t, vector2<double> state, std::string fileName) {
    std::ofstream data;
    data.open(fileName, std::ofstream::app);
    data << t << " " << state.x << " " << state.y << std::endl;
    data.close();
};

// ---------------------------------------------------------------------------------------------------------

// callback to do nothing
auto Nothing = [&](double, vector2<double>, std::string) {};

// ---------------------------------------------------------------------------------------------------------

// main function
int main(int, char **)
{

    // ---------------------------------------------------------------------------------------------------------
    // logistic equation simulation
    /*
    // difference of birth and death for logistic equation
    std::vector<double> rVec{-1., -0.5, 0., 0.5, 1};
    // initial conditions
    std::vector<double> x0Vec{0., 0.2, 0.4, 0.6, 0.8, 1.};
    
    // simulations
    for (int j = 0; j < static_cast<int>(rVec.size()); j++)
    {
        // set parameter
        r = rVec[j];
        for (int i = 0; i < static_cast<int>(x0Vec.size()); i++)
        {
            // set parameter
            double x0 = x0Vec[i];

            // simulation of the logistic equation
            SolverEuler(x0, t0, t1, h, LogisticEq, toFile1D, "../LogE/LogE_r" + std::to_string(j) + "x0" + std::to_string(i) + ".txt");
            SolverRK4(x0, t0, t1, h, LogisticEq, toFile1D, "../LogRK/LogRK_r" + std::to_string(j) + "x0" + std::to_string(i) + ".txt");
        }
    }
    */
    // ---------------------------------------------------------------------------------------------------------
    /*
    // coupled logistic equation simulation
    // initial conditions
    vector2<double> y0{10., 10.};
    // size of square lattice
    int N = 100;
    // container for population numbers
    std::vector<vector2<double>> nContainer(N * N);
    // simulations
    k1 = 100, k2 = 200;
    for (int i = 0; i < N; i++)
    {
        alpha = 0.1 + i * 0.1;
        for (int j = 0; j < N; j++)
        {
            beta = 0.1 + j * 0.1;
            nContainer[j + N * i] = SolverRK4(y0, t0, t1, h, LogisticCoupledEq, Nothing, "None");
        }
    }

    // save one for presentation
    k1 = 100, k2 = 200, alpha = 0.25, beta = 1;
    SolverRK4(y0, t0, t1, h, LogisticCoupledEq, toFile2D, "LogCExample.txt");

    // write results to file
    std::ofstream data;
    data.open("LogCLattice1.txt");
    for (int i = 0; i < N * N; i++)
    {
        data << nContainer[i] << std::endl;
    }
    */
}
