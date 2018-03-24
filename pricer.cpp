#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <iomanip>
#include <Eigen/QR>

using namespace std;
using namespace Eigen;

float phi(float y);
float bs(float values[]);
float fdm(float values[]);

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

    float x;
    float params[11];

    for (int i=0; i<11; i++)
    {
        options >> x;
        params[i] = x;
    }

    // close the parameters file
    options.close();

    float price1 = bs(params);
    float price2 = fdm(params);

    cout << price1 << ", " << price2 << endl;

}

// compute option price using Black-Scholes formulas
float bs(float values[])
{
    // give the variables names that are easier to work with
    int CP = values[2];
    float S = values[3];
    float K = values[4];
    float tau = values[5];
    float sigma = values[6];
    float r = values[7];
    float a = values[8];

    // declare the constants used in the B-S formula
    float d1 = (1/(sigma*sqrt(tau)))*(log(S/K) + (r-a+(pow(sigma,2)/2))*tau);
    float d2 = (1/(sigma*sqrt(tau)))*(log(S/K) + (r-a-(pow(sigma,2)/2))*tau);
    float price;

    if (CP == 0) {
        price = S*exp(-a*tau)*phi(d1) - K*exp(-r*tau)*phi(d2);
    }

    else {
        price = K*exp(-r*tau)*phi(-d2) - S*exp(-a*tau)*phi(-d1);
    }

    return price;
}

// helper functions for fdm
float aj(int j, int M, float sig, float r, float q, float t) {
    return -0.5*((r-q)*(M-j)*t + pow(sig*(M-j), 2)*t);
}

float bj(int j, int M, float sig, float r, float q, float t) {
    return 1 + pow(sig*(M-j),2)*t + r*t;
}

float cj(int j, int M, float sig, float r, float q, float t) {
    return 0.5*((r-q)*(M-j)*t - pow(sig*(M-j),2)*t);
}

// finite difference method
float fdm(float values[])
{
    // be careful when using M and N to make sure that floats are returned
    // when appropriate
    float EA = values[1];
    float CP = values[2];
    float S = values[3];
    float K = values[4];
    float T = values[5];
    float sig = values[6];
    float r = values[7];
    float q = values[8];
    int M = values[9];
    int N = values[10];

    // mult will be used to determine the maximum stock price of the grid
    int mult = 5;

    // change M slightly to make it easier to recover option price at the end
    M -= (M % mult);
    float Smax = mult*S;

    float s = Smax/M;
    float t = T/N;

    // initialize the grid with all zeros
    MatrixXf g(M+1, N+1);
    g.setZero();

    // set the boundary conditions corresponding to whether the option
    // is a call or put
    for (int i=0; i<N+1; i++) {
        g(0, i) += (1-CP)*(Smax - K*exp(-(r-q)*(T-t*i)));
        g(M, i) += CP*K*exp(-(r-q)*(T-t*i));
    }

    for (int i=1; i<M; i++) {
        g(i, N) = (1-CP)*std::max((M-i)*s-K, (float)0.0) + CP*std::max(K-(M-i)*s, (float)0.0);
    }

    // make sure the row corresponding to the current price actually
    // contains the current price in the first column and the appropriate
    // intrinsic value in the last column
    g(M - M/mult, 0) = S;
    g(M - M/mult, N) = (1-CP)*std::max(S-K, (float)0.0) + CP*std::max(K-S, (float)0.0);

    /* Here we update the grid column-by-column, moving right to left. We have
    to solve a system of linear equations for each column, but we can perform
    back substitution on the coefficient matrices to obtain the solutions,
    which is much faster than inverting them. Think of A below as follows:
    A = [B, b], where B is the square coefficient matrix and b is the vector
    of right-hand side values. We are solving to obtain x = B^{-1}b. */

    MatrixXf A(M+1, M+1);
    VectorXf d(M+1);
    VectorXf x(M+1);

    for (int i=N; i>0; i--) {

        // set A equal to identity
        A = MatrixXf::Identity(M+1,M+1);
        d = g.col(i);

        for (int j=1; j<M; j++) {
            A(j,j-1) = aj(j, M, sig, r, q, t);
            A(j,j) = bj(j, M, sig, r, q, t);
            A(j,j+1) = cj(j, M, sig, r, q, t);
        }

        // fast solution of equation Ax=d using special structure of A
        double a;

        // zero out the first column, which has only one other non-zero entry,
        // and update RHS accordingly
        d(1) -= A(1,0)*d(0);
        A(1,0) = 0;

        // normalize i,i-th entry, zero out everything below it, update RHS
        for (int j=1; j<M; j++) {
            a = 1/A(j,j);
            A(j,j) *= a;
            A(j,j+1) *= a;
            d(j) *= a;

            A(j+1,j+1) -= A(j+1,j)*A(j,j+1);
            d(j+1) -= A(j+1,j)*d(j);

            A(j+1,j) = 0;
        }

        // start process over, this time zeroing out the entries above the diagonal
        d(M-1) -= A(M-1,M)*d(M);
        A(M-1,M) = 0;

        // to obtain the correct b we don't need to set A(i-1,i)=0, but since
        // we are reusing A we might as well reset it to the identity
        for (int j=M-1; j>0; j--) {
            d(j-1) -= A(j-1,j)*d(j);
            A(j-1,j) = 0;
        }

        for (int j=0; j<M+1; j++) {
            d(j) = (1-CP)*std::max(d(j), EA*((M-j)*s - K)) + CP*(std::max(d(j), EA*(K-(M-j)*s)));
        }

        // print out quantities in Ax = d to use in debugging
        // cout << endl << A << endl << d << endl << x << endl;

        g.col(i-1) = d;

        // these two lines are inefficient and unnecessary -- find out how
        // to get rid of them
        g(0, i-1) = (1-CP)*(Smax - K*exp(-(r-q)*(T-t*i)));
        g(M, i-1) = CP*K*exp(-(r-q)*(T-t*i));
    }

    return g(M - M/mult, 0);
}

// copied this directly from John Cook at https://www.johndcook.com/blog/cpp_phi/
float phi(float x)
{
    // constants
    float a1 =  0.254829592;
    float a2 = -0.284496736;
    float a3 =  1.421413741;
    float a4 = -1.453152027;
    float a5 =  1.061405429;
    float p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    float t = 1.0/(1.0 + p*x);
    float y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}
