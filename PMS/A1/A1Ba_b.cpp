#include <iostream>
#include <cmath>
typedef long long ll;
using namespace std;
// Define the function f(x)
double f(double x, double k, double z)
{
  return (k - 1) * z / (x * (k - 1) + 1);
}

// Define the derivative f'(x)
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
  return (k)*z / (x * (k - 1) + 1);
}
int main()
{
  {
    float F = 100;
    float z[3];
    z[0] = 0.6, z[1] = 0.25, z[2] = 0.15;

    float Pt = 1.01325; // converted to bar;
    cout << "Input temperature in Celcius" << endl;
    float temp;
    cin >> temp;
    temp += 273.18;
    // pressure and equilibrium constant
    float P[3], K[3];
    // index 0 is benzene
    // index 1 is toluene
    // index 2 is o-xylene
    float A[3] = {4.01814, 4.14157, 4.12928};
    float B[3] = {1203.835, 1377.578, 1478.244};
    float C[3] = {-53.226, -50.507, -59.076};

    for (int i = 0; i < 3; i++)
    {
      float t = A[i] - (B[i] / (C[i] + temp));
      P[i] = pow(static_cast<double>(10), t);

      K[i] = P[i] / Pt;
    }
    for (int i = 0; i < 3; i++)
    {
      cout << K[i] << endl;
    }

    double x, prev;

    int id = 0;
    double tol = 0.000001;
    bool flag = true;
    double a = 0.0, b = 1.0, fa = 0, fb = 0;
    for (int i = 0; i < 3; i++)
    {
      fa += f(a, K[i], z[i]);
      fb += f(b, K[i], z[i]);
    }

    if (fa * fb > 0)
    {
      cout << "No valid root found in [0,1]. Check K-values and change temperature." << endl;
      return -1;
    }
    cout << "Choose initial guess:" << endl;
    cin >> x;
    while (id < 100)
    {
      float summation = 0, summationed = 0;
      for (int i = 0; i < 3; i++)
      {
        summation += f(x, K[i], z[i]);
        summationed += df(x, K[i], z[i]);
        // x -= fx / dfx
        if (id % 5 == 0)
          cout << f(x, K[i], z[i]) << " " << df(x, K[i], z[i]) << endl;
      }
      if (summationed == 0)
      {
        cout << id << endl;
        cout << "Choose other initial guess.Try again" << endl;
        flag = false;
        break;
      }

      if (tol >= abs(x - prev))
      {
        cout << "value of x in " << id << " iteration is: " << x << endl;
        break;
      }
      if (id + 1 >= 100)
      {
        flag = false;
        cout << "Max iterations reached. Solution may not be accurate.  " << x << endl;
      }
      prev = x;
      x -= summation / summationed;
      id++;
    }
    if (flag)
    {

      float V = x * F;
      float L = F - V;
      float X[3], Y[3];
      for (int i = 0; i < 3; i++)
      {
        X[i] = ddx(x, K[i], z[i]);
        Y[i] = ddy(x, K[i], z[i]);
      }
      cout << "Vapour flow rate is: " << V << endl;
      cout << "Liquid flow rate is: " << L << endl;
      cout << "The component of Benzene is liquid phase and vapour phase respectively are: " << endl;
      cout << X[0] << " " << Y[0] << endl;
      cout << "The component of Toulene is liquid phase and vapour phase respectively are: " << endl;
      cout << X[1] << " " << Y[1] << endl;
      cout << "The component of o-Xylene is liquid phase and vapour phase respectively are: " << endl;
      cout << X[2] << " " << Y[2] << endl;
    }
  }
}
