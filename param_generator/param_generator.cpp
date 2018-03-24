#include <iostream>
#include <fstream>
#include "lcgrand.h"

using namespace std;

int main() {
    // want to generate option parameters to feed into pricer.cpp

    /* 0) 0 for B-S, 1 for FDM
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

    // open grid dimensions textfile and copy it into an array
    ifstream grid_dims;
    grid_dims.open("grid_dimensions.txt");

    if (!grid_dims) {
        std::cerr << "Unable to open grid_dimensions.txt. \n";
        exit(1);
    }

    int dimension[50];
    for (int k=0; k<50; k++) {
        grid_dims >> dimension[k];
    }

    ofstream options;
    options.open("option_params.txt");

    if (!options) {
        std::cerr << "Unable to open option_params.txt. \n";
        exit(1);
    }

    float S, K, T, sig, r, q;
    int M, N, seed = 1;
    int ops = 20; // number of options to generate

    for (int i=0; i<ops; i++) {
        S = 20 + 400*lcgrand(seed); // generate current underlying price
        K = S + min((float)100, S)*(-1+2*lcgrand(seed)); // strike price within min of
            // 100 and S of S
        T = lcgrand(seed); // expiry sometime in the next year
        sig = 0.1 + 0.9*lcgrand(seed); // volatility between 0.1 and 1
        r = 0.01 + 0.19*lcgrand(seed); // risk-free rate between 0.01 and 0.2
        q = r*lcgrand(seed); // dividend rate between 0 and r

        for (int k=0; k<50; k++) {
            M = dimension[k]; N = dimension[k];
            for (int j=0; j<2; j++) {
                options << 0 << " " // BS, FDM designator; not used at the moment
                    << 0 << " " // European
                    << j << " " // 0 for call, 1 for put
                    << S << " " // current price
                    << K << " " // strike
                    << T << " " // expiry
                    << sig << " " // volatility
                    << r << " " // risk-free rate
                    << q << " " // dividend rate
                    << M << " " // number of price steps
                    << N; // number of time steps

                    if (i==ops-1 & j==1 & k==49) {
                        options << "*";
                    }
                    else {
                        options << "\n";
                    }
            }
        }
    }

    return 0;
}
