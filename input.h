//
// Created by epsilon on 03.08.24.
//

#ifndef COMPMOD_INPUT_H
#define COMPMOD_INPUT_H

#include <complex>

const long double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

std::string outExactFilename = "exact.csv";
std::string outCalcFilename = "calculated.csv";
std::string outDiffFilename = "difference.csv";


long double f_a(long double t, long double x, long double u) {
    return u;
}

long double df_a(long double t, long double x, long double u) {
    return 1;
}

long double f_beg(long double x) {
    return 0.25 + 0.5*std::sin(PI*x);
}

long double f_beg_disc(long double x) {
    return x >= 0l ? 1 : -1;
}

long double X1 = -1;
long double X2 = 1;
long double T1 = 0;
long double T2 = 2;
long double CFL = 0.9l;
unsigned int n_x = 100;


/*
long double f_hoppf(long double t, long double x) {
    return 1;
}

long double f_beg2(long double x) {
    return std::exp(-5.0l*(x-0.5)*(x-0.5));
}

long double f_beg_disc(long double x) {
    return x >= 0.5l && x < 1 ? 1 : -1;
}

long double f_beg_disc2(long double x) {
    return x >= 0.25l && x < 0.75 ? 1 : 0;
}
*/
#endif //COMPMOD_INPUT_H