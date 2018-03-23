#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <iomanip>
#include <algorithm>
#include <Eigen/LU>

using namespace std;
using namespace Eigen;

double phi(double y);
double bs(double values[]);
double fdm(double values[]);

int main()
{
    // read in option parameters file
    std::ifstream options;
    options.open("option_params.txt");

    // check if file containing option parameters was read in successfully
    if (!options) {
        std::cerr << "Unable to open option_params.txt. \n";
        exit(1);
    }

    /* read first seven elements of option_params into a vector
    parameters are, in order:
    0) 0 for B-S, 1 for FDM
    1) 0 for European, 1 for American;
    2) 0 for call, 1 for put;
    3) current price;
    4) strike;
    5) time until expiration;
    6) volatility;
    7) risk-free interest rate;
    8) continuous dividend rate
    9) number of price steps
    10) number of time steps */

    double x;
    double params[11];

    for (int i=0; i<11; i++)
    {
        options >> x;
        params[i] = x;
    }

    // close the parameters file
    options.close();

    double price1 = bs(params);
    double price2 = fdm(params);

    cout << price1 << ", " << price2 << "\n";
}

// compute option price using Black-Scholes formulas
double bs(double values[])
{
    // give the variables names that are easier to work with
    int CP = values[2];
    double S = values[3];
    double K = values[4];
    double tau = values[5];
    double sigma = values[6];
    double r = values[7];
    double a = values[8];

    // declare the constants used in the B-S formula
    double d1 = (1/(sigma*sqrt(tau)))*(log(S/K) + (r-a+(pow(sigma,2)/2))*tau);
    double d2 = (1/(sigma*sqrt(tau)))*(log(S/K) + (r-a-(pow(sigma,2)/2))*tau);
    double price;

    if (CP == 0) {
        price = S*exp(-a*tau)*phi(d1) - K*exp(-r*tau)*phi(d2);
    }

    else {
        price = K*exp(-r*tau)*phi(-d2) - S*exp(-a*tau)*phi(-d1);
    }

    return price;
}

// helper functions for fdm
double aj(int j, int M, float sig, float r, float q, float t) {
    return -0.5*((r-q)*(M-j)*t - pow(sig*(M-j), 2)*t);
}

double bj(int j, int M, float sig, float r, float q, float t) {
    return 1 + pow(sig*(M-j),2)*t + r*t;
}

double cj(int j, int M, float sig, float r, float q, float t) {
    return 0.5*((r-q)*(M-j)*t - pow(sig*(M-j),2)*t);
}

// finite difference method
double fdm(double values[])
{
    // be careful when using M and N to make sure that doubles are returned
    // when appropriate
    double EA = values[1];
    double CP = values[2];
    double S = values[3];
    double K = values[4];
    double T = values[5];
    double sig = values[6];
    double r = values[7];
    double q = values[8];
    int M = values[9];
    int N = values[10];

    // mult will be used to determine the maximum stock price of the grid
    int mult = 5;

    // change M slightly to make it easier to recover option price at the end
    M -= (M % mult);
    double Smax = mult*S;

    double s = Smax/M;
    double t = T/N;

    // initialize the grid with all zeros
    double g[M+1][N+1];
    for (int i=0; i<M+1; i++) {
        for (int j=0; j<N+1; j++) {
            g[i][j] = 0;
        }
    }

    // set the boundary conditions corresponding to whether the option
    // is a call or put
    for (int i=0; i<N+1; i++) {
        g[0][i] += (1-CP)*(Smax - K*exp(-(r-q)*(T-t*i)));
        g[M][i] += CP*K*exp(-(r-q)*(T-t*i));
    }

    for (int i=1; i<M; i++) {
        g[i][N] = (1-CP)*std::max((M-i)*s-K, 0.0) + CP*std::max(K-(M-i)*s, 0.0);
    }

    // make sure the row corresponding to the current price actually
    // contains the current price in the first column and the appropriate
    // intrinsic value in the last column
    g[M - M/mult][0] = S;
    g[M - M/mult][N] = (1-CP)*std::max(S-K, 0.0) + CP*std::max(K-S, 0.0);

    /* Here we update the grid column-by-column, moving right to left. We have
    to solve a system of linear equations for each column, but we can perform
    back substitution on the coefficient matrices to obtain the solutions,
    which is much faster than inverting them. Think of A below as follows:
    A = [B, b], where B is the square coefficient matrix and b is the vector
    of right-hand side values. We are solving to obtain x = B^{-1}b. */

    MatrixXd A(M+1, M+1);
    VectorXd d(M+1);
    VectorXd x(M+1);

    for (int i=N; i>0; i--) {

        // set A equal to identity
        A = MatrixXd::Identity(M+1,M+1);

        for (int j=0; j<M+1; j++) {
            d(j) = g[j][i];
        }

        for (int j=1; j<M; j++) {
            A(j,j-1) = aj(j, M, sig, r, q, t);
            A(j,j) = bj(j, M, sig, r, q, t);
            A(j,j+1) = cj(j, M, sig, r, q, t);
        }

        x = A.partialPivLu().solve(d);

        for (int j=0; j<M+1; j++) {
            x(j) = (1-CP)*std::max(x(j), EA*((M-j)*s - K)) + CP*(std::max(x(j), EA*(K-(M-j)*s)));
        }

        for (int j=1; j<M; j++) {
            g[j][i-1] = x(j);
        }

    }

    return g[M - M/mult][0];
}

// copied this directly from John Cook at https://www.johndcook.com/blog/cpp_phi/
double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}
