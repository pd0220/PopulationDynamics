// population dynamics simulation implementations

// used headers and libraries
#include <iostream>
#include "NumericalIntegration.hh"
#include "Vector2D.hh"

// difference of birth and death for logistic equation
double const r = 0.5;

// differences of birth and death for coupled logistic equations
double const r1 = 0.5;
double const r2 = 0.5;
// coupling coefficients for coupled logistic equations
double const alpha = 1.;
double const beta = 1.;
// saturation values for coupled logistic equations
double const k1 = 1.;
double const k2 = 1.;

// parameters for the Lotka-Volterra model
double const a = 20.;
double const b = 1.;
double const c = 30.;
double const d = 1.;

// RHS for logistic equation
auto LogisticEq = [&](double, double x) {
    return r * x * (1 - x);
};

// RHS for coupled logistic equations
auto LogisticCoupledEq = [&](double, vector2<double> state) -> vector2<double> {
    return {
        r1 * state.x * (1 - (state.x + alpha * state.y) / k1),
        r2 * state.y * (1 - (state.y + beta * state.x) / k2)};
};

// RHS for Lotka-Volterra model
auto LotkaVolterraEq = [&](double, vector2<double> state) -> vector2<double> {
    return {
        a * state.x - b * state.x * state.y,
        c * state.x * state.y - d * state.y};
};

// main function
int main(int, char **)
{
}
