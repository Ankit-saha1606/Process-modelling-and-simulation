#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

// Function definitions
double f(double x, double k, double z)
{
    return (k - 1) * z / (x * (k - 1) + 1);
}

double df(double x, double k, double z)
{
    return -(k - 1) * (k - 1) * z / ((x * (k - 1) + 1) * (x * (k - 1) + 1));
}

double ddx(double x, double k, double z)
{
    return z / (x * (k - 1) + 1);
}

double ddy(double x, double k, double z)
{
    return k * z / (x * (k - 1) + 1);
}

// Rachford-Rice using Newton-Raphson
double solveForX(double K[], double Z[], int components, double tol = 1e-6)
{
    double x = 1; // Initial guess
    int max_iter = 100;

    for (int iter = 0; iter < max_iter; iter++)
    {
        double sum_f = 0, sum_df = 0;
        for (int i = 0; i < components; i++)
        {
            sum_f += f(x, K[i], Z[i]);
            sum_df += df(x, K[i], Z[i]);
        }
        // cout << x << " " << sum_f << endl;
        if (fabs(sum_df) < tol)
        {
            cout << "Derivative too small, try different initial guess." << endl;
            return -1;
        }

        double x_new = x - sum_f / sum_df;
        if (fabs(x_new - x) < tol)
        {
            return x_new;
        }

        x = x_new;
    }

    cout << "Newton-Raphson did not converge." << endl;
    return -1;
}

// Flash Distillation Calculation
void flash_distillation(double F, double temp, double &V, double &L, double X[], double Y[], double z[])
{
    double Pt = 1.01325; // Pressure in bar

    // Antoine coefficients
    double A[3] = {4.01814, 4.14157, 4.12928};
    double B[3] = {1203.835, 1377.578, 1478.244};
    double C[3] = {-53.226, -50.507, -59.076};

    double P[3], K[3];

    // Calculate K-values
    for (int i = 0; i < 3; i++)
    {
        double t = A[i] - (B[i] / (C[i] + temp + 273.18));
        P[i] = pow(10, t);
        K[i] = P[i] / Pt;
    }

    // Solve Rachford-Rice Equation
    // cout << z[0] << endl;
    double x = solveForX(K, z, 3);
    if (x == -1)
    {
        cout << "Rachford-Rice method failed. Adjust parameters." << endl;
        return;
    }

    // Compute Vapor & Liquid Flow Rates
    V = x * F;
    L = F - V;

    // Compute phase compositions
    for (int i = 0; i < 3; i++)
    {
        X[i] = ddx(x, K[i], z[i]); // Liquid phase
        Y[i] = ddy(x, K[i], z[i]); // Vapor phase
    }
}

int main()
{
    double F = 100, temp; // Initial feed & temperature
    cout << "Give initial temperature between 91 and 108 degree Celsius: ";
    cin >> temp;

    vector<double> benzene_recovery; // Store benzene recovery for each cycle

    double V, L, X[3], Y[3];
    double Z[3] = {0.6, 0.25, 0.15}; // Initial feed composition

    for (int cycle = 0; cycle < 10; cycle++)
    {
        flash_distillation(F, temp, V, L, X, Y, Z);

        double benzene_in_vapor = Y[0] * V; // Recovered benzene
        benzene_recovery.push_back(benzene_in_vapor);

        cout << "Cycle " << cycle + 1 << ": Benzene in Vapor = " << benzene_in_vapor << " kmol" << endl;

        // Stop if all liquid is vaporized
        if (L < 1e-6)
        {
            cout << "Stopping at cycle " << cycle + 1 << " as no liquid remains.\n";
            break;
        }

        // Update new feed as the recycled bottom product
        F = L;

        // Update the new feed mole composition as the actual moles (not just fractions)
        Z[0] = X[0] * L;
        Z[1] = X[1] * L;
        Z[2] = X[2] * L;

        // Normalize Z to prevent mole fraction drift
        double total = Z[0] + Z[1] + Z[2];
        Z[0] /= total;
        Z[1] /= total;
        Z[2] /= total;
    }
    double sumbenzene = 0;
    for (int i = 0; i < 10; i++)
    {
        sumbenzene += benzene_recovery[i];
    }
    // Print final benzene recovery after 100 cycles
    cout << "\nTotal Benzene Recovered after " << benzene_recovery.size() << " cycles: "
         << sumbenzene << " kmol" << endl;

    return 0;
}
