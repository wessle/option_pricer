#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <iomanip>

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

    double price = bs(params);

    std::cout << std::setprecision(10) << price << "\n";
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

double fdm(double values[])
{
    
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
