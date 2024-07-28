#include <iostream>
#include <cmath>
#include "solver_fvm.h"
#include "verifier.h"

long double f_a(long double t, long double x) {
    return 1;
}

long double f_beg(long double x) {
    return sin(2*3.1415*x);
}


void print_vector2d(const std::vector<std::vector<long double>> &field) {
    for (int n = 0; n < field.size(); n++) {
        for (int i = 0; i < field[0].size(); i++) {
            printf("%Lf,", field[n][i]);
        }
        printf("\n");
    }
}

int main() {
    long double T, tau, X, h;
    T = 2;
    tau = 0.01;
    X = 1;
    h = 0.01;

    SolverGodunov s1 = SolverGodunov(T, tau, X, h, &f_a, &f_beg);
    s1.calculate();

    std::vector<std::vector<long double>> field = s1.getField();

    Verifier v1 = Verifier(field, T, tau, X, h, f_a(0, 1), &f_beg);
    print_vector2d(v1.getDiff());
    printf("%Lf", v1.get_L1());


}
