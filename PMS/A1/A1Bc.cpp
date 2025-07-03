#include <iostream>
#include <cmath>
#include <limits>

typedef long long ll;
using namespace std;

// Define the function f(x) for Rachford-Rice equation
double f(double x, double k, double z)
{
    return (k - 1) * z / (x * (k - 1) + 1);
}

// Define the derivative f'(x)
double df(double x, double k, double z)
{
    return -(k - 1) * (k - 1) * z / pow((x * (k - 1) + 1), 2);
}

// Compute liquid phase composition
double ddx(double x, double k, double z)
{
    return z / (x * (k - 1) + 1);
}

// Compute vapor phase composition
double ddy(double x, double k, double z)
{
    return (k)*z / (x * (k - 1) + 1);
}

// Solve for x using Newton-Raphson method
double solveForX(double K[], double Z[], int components, double tol = 1e-6)
{
    double x = 0.5; // Initial guess
    int max_iter = 100;
    double a = 0.0, b = 1.0, fa = 0, fb = 0;
    for (int i = 0; i < 3; i++)
    {
        fa += f(a, K[i], Z[i]);
        fb += f(b, K[i], Z[i]);
    }

    if (fa * fb > 0)
    {
        // cout << "No valid root found in [0,1]. Check K-values and change temperature." << endl;
        return -1;
    }
    for (int iter = 0; iter < max_iter; iter++)
    {
        double sum_f = 0, sum_df = 0;
        for (int i = 0; i < components; i++)
        {
            sum_f += f(x, K[i], Z[i]);
            sum_df += df(x, K[i], Z[i]);
        }

        if (fabs(sum_df) < tol)
        { // Avoid division by zero
          // cout << "Derivative too small, try different initial guess." << endl;
            return -1;
        }

        double x_new = x - sum_f / sum_df;

        if (fabs(x_new - x) < tol)
        { // Convergence check
            return x_new;
        }

        x = x_new;
    }

    cout << "Newton-Raphson did not converge." << endl;
    return -1;
}

int main()
{
    float F = 100;                   // Feed flow rate
    double Z[3] = {0.6, 0.25, 0.15}; // Initial mole fractions

    float Pt = 1.01325; // Pressure in bar

    // Antoine constants for Benzene, Toluene, o-Xylene
    float A[3] = {4.01814, 4.14157, 4.12928};
    float B[3] = {1203.835, 1377.578, 1478.244};
    float C[3] = {-53.226, -50.507, -59.076};

    double bestTemp = 0, maxRecovery = -1; // Store optimal values
    // it can be shown in previous question that the temperature limit is 108 degree Celcius
    for (float tempC = 100; tempC <= 108; tempC += 0.0001)
    {
        double temp = tempC + 273.18; // Convert to Kelvin

        double P[3], K[3];

        // Compute equilibrium constants K
        for (int i = 0; i < 3; i++)
        {
            float t = A[i] - (B[i] / (C[i] + temp));
            P[i] = pow(10, t);
            K[i] = P[i] / Pt;
        }

        // Solve for vapor fraction x using Newton-Raphson
        double x = solveForX(K, Z, 3);

        if (x == -1)
            continue; // Skip if solution is invalid

        double V = x * F;
        double L = F - V;
        double X[3], Y[3];

        for (int i = 0; i < 3; i++)
        {
            X[i] = ddx(x, K[i], Z[i]);
            Y[i] = ddy(x, K[i], Z[i]);
        }

        // Compute benzene recovery
        double recoveryBenzene = (V * Y[0]) / (F * Z[0]);

        // Update best temperature
        if (recoveryBenzene > maxRecovery)
        {
            maxRecovery = recoveryBenzene;
            bestTemp = tempC;
        }
    }

    cout << "Optimal temperature for maximum benzene recovery: " << bestTemp << "°C" << endl;
    cout << "Maximum benzene recovery: " << maxRecovery * 100 << "%" << endl;

    return 0;
}
