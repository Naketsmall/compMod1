//
// Created by epsilon on 28.07.24.
//

#include "solver_fvm.h"


long double SolverFVM::get_courant(double t, double x) {
    return this->f_a(t, x) * this->tau / this->h;
}

void SolverFVM::calculate_cell(int n, int i) {
    field[n][i] = field[n - 1][i] - tau / h * (w_r(n - 1, i) - w_l(n - 1, i));
}

SolverFVM::SolverFVM(long double T, long double tau, long double X, long double h,
                     long double (*f_a)(long double, long double),
                     long double (*f_beg)(long double)) : T(T), tau(tau), X(X), h(h), f_a(f_a) {

    n_t = (int) (T / tau)+1;
    n_x = (int) (X / h)+1;

    for (int n = 0; n < n_t; n++) {
        field.emplace_back(n_x);
    }
    for (int i = 0; i < n_x; i++) {
        field[0][i] = f_beg(i*h);
    }

}

int SolverFVM::calculate() {
    for (int n = 1; n < n_t; n++) {
        for (int i = 0; i < n_x; i++) {
            calculate_cell(n, i);
        }
    }
}

const std::vector<std::vector<long double>> &SolverFVM::getField() const {
    return field;
}



long double SolverGodunov::w_l(int n, int i) {
    return (f_a(n * tau, i * h) + f_a(n * tau, (i - 1) * h)) / 2 * field[n][(n_x + i - 1) % n_x];
}

long double SolverGodunov::w_r(int n, int i) {
    return (f_a(n * tau, i * h) + f_a(n * tau, (i + 1) * h)) / 2 * field[n][i];
}

SolverGodunov::SolverGodunov(long double T, long double tau, long double X, long double h,
                             long double (*f_a)(long double, long double), long double (*f_beg)(long double))
        : SolverFVM(T, tau, X, h, f_a, f_beg) {

}
