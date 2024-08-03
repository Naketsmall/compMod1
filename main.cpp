#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "solver_fvm.h"
#include "verifier.h"

/*
 *
 */

const long double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

long double f_a(long double t, long double x, long double u) {
    return 1;
}

long double f_hoppf(long double t, long double x) {
    return 1;
}

long double f_beg(long double x) {
    return std::sin(2*PI*x);
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


std::string vector2d_to_string(const std::vector<std::vector<long double>> &field) {
    std::stringstream line;

    for (int n = 0; n < field.size(); n++) {
        for (int i = 0; i < field[0].size(); i++) {
            if (i != field[0].size()-1)
                line << std::fixed << std::setprecision(30) << field[n][i] << ",";// field[n][i] + ",";
            else
                line << std::fixed << std::to_string(field[n][i]);
        }
        line << "\n";
    }
    return std::move(line).str();
}

int main() {
    long double T, tau, X, h, CFL;
    X = 1;
    h = 1.0l/3200;
    T = 2;
    CFL = 0.8l;

    tau = CFL * h * f_a(0,0, 0);

    SolverGodunov s1 = SolverGodunov(T, tau, X, h, &f_a, &f_beg_disc2);
    s1.calculate();

    std::vector<std::vector<long double>> field = s1.getField();

    Verifier v1 = Verifier(field, T, tau, X, h, f_a(0, 1, 0), &f_beg_disc2);


    std::ofstream out_acc("accurate11.csv");
    out_acc << vector2d_to_string(v1.getAccurate());
    out_acc.close();

    std::ofstream out_calc("calculated11.csv");
    out_calc << vector2d_to_string(field);
    out_calc.close();

    std::ofstream out_diff("diff11.csv");
    out_diff << vector2d_to_string(v1.getDiff());
    out_diff.close();

    printf("Max: %.30Lf\n", v1.get_max());
    printf("L1: %.30Lf\n", v1.get_L1());
    printf("shape calc: %zu, %zu\n", field.size(), field[0].size());
    printf("shape acc: %zu, %zu\n", v1.getAccurate().size(), v1.getAccurate()[0].size());
    printf("CFL: %.30Lf\n", s1.get_courant(0, 0, 0));


}
