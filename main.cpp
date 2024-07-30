#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "solver_fvm.h"
#include "verifier.h"

const long double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

long double f_a(long double t, long double x) {
    return 1;
}

long double f_beg(long double x) {
    return std::sin(2*PI*x);
}

long double f_beg2(long double x) {
    return std::exp(-5.0l*(x-0.5)*(x-0.5));
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
    long double T, tau, X, h;
    T = 2;
    tau = 0.01;
    X = 1;
    h = 0.01;

    SolverGodunov s1 = SolverGodunov(T, tau, X, h, &f_a, &f_beg2);
    s1.calculate();

    std::vector<std::vector<long double>> field = s1.getField();

    Verifier v1 = Verifier(field, T, tau, X, h, f_a(0, 1), &f_beg2);


    std::ofstream out_acc("accurate2.csv");
    out_acc << vector2d_to_string(v1.getAccurate());
    out_acc.close();

    std::ofstream out_calc("calculated2.csv");
    out_calc << vector2d_to_string(field);
    out_calc.close();

    std::ofstream out_diff("diff2.csv");
    out_diff << vector2d_to_string(v1.getDiff());
    out_diff.close();

    printf("L1: %.30Lf\n", v1.get_L1());
    printf("Max: %.30Lf\n", v1.get_max());
    printf("shape calc: %zu, %zu\n", field.size(), field[0].size());
    printf("shape acc: %zu, %zu\n", v1.getAccurate().size(), v1.getAccurate()[0].size());


}
