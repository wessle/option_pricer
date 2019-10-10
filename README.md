# option_pricer
Option pricing software project. The main option pricing program is written in C++ and able to price European and American options using the Black-Scholes formula and the method of finite differences. Python files are prototypes, but are more user friendly -- use them in the interactive interpreter to calculate prices from the command line.

The subdirectories in this directory contain various versions of the pricing engine. Each version reads in a set of option parameters from a text file, performs a series of computations, then outputs to a separate text file. I used these programs to perform numerical experiments to determine the effectiveness of my implementation of the implicit finite difference method for calculating option prices, then output the resulting data in a form that is easy to translate into graphs and other figures.

To take a peak at the main elements of the pricing engine, check out pricer.cpp. pricer.cpp reads in sets of option parameters from 'option_params.txt', computes the Black-Scholes and finite difference prices, then outputs the option parameters and error between the two prices to 'errors.txt'. The functions 'bs()' and 'fdm()' are declared and used in pricer.cpp. 'bs()' contains my Black-Scholes implementation, while 'fdm()' contains my finite difference method implementation. 'bs()' and 'fdm()' are the core of this project.

If you want to compile pricer.cpp, you must compile it with pricer.h. Most C++ files in the subdirectories are versions of pricer.cpp -- they generally must be compiled with their corresponding pricer.h versions.

All C++ files were compiled with C++11 as standard.
