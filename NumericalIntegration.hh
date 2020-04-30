// numerical solvers for differential equations

// used headers and libraries
#include <fstream>

// simple Euler-method
template <typename State, typename T, typename RHS, typename Callback>
auto SolverEuler(State y0, T t0, T t1, T h, RHS f, Callback cb)
{
    // initial conditions
    T t = t0;
    State y = y0;
    // integration
    while (t < t1)
    {
        // do not overflow in time
        if (t + h > t1)
            h = t1 - t;
        // Euler step
        y += h * f(t, y);
        t += h;
        // callback function
        cb(t, y);
    }
}

// 4th order Runge-Kutta method
template <typename State, typename T, typename RHS, typename Callback>
auto SolverRK4(State y0, T t0, T t1, T h, RHS f, Callback cb)
{
    // initial conditions
    T t = t0;
    State y = y0;
    // integration
    while (t < t1)
    {
        // do not overflow in time
        if (t + h > t1)
            h = t1 - t;
        // RK4 step
        State k1 = f(t, y);
        State k2 = f(t + h * (T)0.5, y + (h * (T)0.5) * k1);
        State k3 = f(t + h * (T)0.5, y + (h * (T)0.5) * k2);
        State k4 = f(t + h, y + h * k3);
        y += (k1 + k4 + (T)2 * (k2 + k3)) * h / (T)6;
        t += h;
        // callback function
        cb(t, y);
    }
}