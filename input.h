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
    return u*u/2;
}

long double df_a(long double t, long double x, long double u) {
    return u;
}

long double f_beg(long double x) {
    return 0.25 + 0.5*std::sin(PI*x);
}

long double f_beg_disc(long double x) {
    return x >= 0l ? 1 : -1;
}

class SolverGodunov : public SolverFVM {
    long double w_r(int n, int i, long double t) override {
        long double u_r;
        if (field[n][i] >= field[n][(n_x + i+1) % n_x]) {
            if (field[n][i] + field[n][(n_x + i+1) % n_x] >= 0)
                u_r = field[n][i];
            else
                u_r = field[n][(n_x + i+1) % n_x];
        } else {
            if (i*h/t <= field[n][i])
                u_r = field[n][i];
            else if (field[n][i] < i*h/t <= field[n][(n_x + i+1) % n_x])
                u_r = i*h/t;
            else
                u_r = field[n][(n_x + i+1) % n_x];
        }
        return f_a(t, h*i, u_r);
        /*
        if (df_a(taus[n], i * h, field[n][i]) + df_a(taus[n], (i + 1) * h, field[n][(i+1) % n_x]) >= 0)
            return (df_a(taus[n], i * h, field[n][i]) + df_a(taus[n], (i + 1) * h, field[n][(i+1) % n_x])) / 2 *
                   field[n][i];
        else
            return (df_a(taus[n], i * h, field[n][i]) + df_a(taus[n], (i + 1) * h, field[n][(i+1) % n_x])) / 2 *
                   field[n][(n_x + i + 1) % n_x];*/
    }

public:
    SolverGodunov(long double T1, long double T2, long double X1, long double X2, unsigned int n_x,
                  long double CFL,
                  long double (*f_a)(long double, long double, long double),
                  long double (*df_a)(long double, long double, long double),
                  long double (*f_beg)(long double)) : SolverFVM(T1, T2, X1, X2, n_x, CFL, f_a, df_a, f_beg) {}


};



long double X1 = -1;
long double X2 = 1;
long double T1 = 0;
long double T2 = 1;
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
