//
// Created by epsilon on 28.07.24.
//

#include <cmath>
#include "solver_fvm.h"


long double SolverFVM::get_courant(long double t, long double x, long double u) {
    return this->f_a(t, x, u) * this->tau / this->h;
}

void SolverFVM::calculate_cell(int n, int i) {
    field[n][i] = field[n - 1][i] - tau / h * (w_r(n - 1, i) - w_l(n - 1, i));
}

SolverFVM::SolverFVM(long double T, long double tau, long double X, long double h,
                     long double (*f_a)(long double, long double, long double),
                     long double (*f_beg)(long double)) : T(T), tau(tau), X(X), h(h), f_a(f_a) {

    n_t = std::ceil(T / tau);
    n_x = std::ceil(X / h);

    for (int n = 0; n < n_t; n++) {
        field.emplace_back(n_x);
    }

    for (int i = 0; i < n_x; i++) {
        field[0][i] = f_beg((i+0.5l)*h);
    }

}
// Необходимо доделать определение через n_x и CFL, чтобы tau рассчитывалось динамически.

SolverFVM::SolverFVM(long double T, int n_t, long double X, int n_x, long double (*f_a)(long double, long double, long double),
                     long double (*f_beg)(long double)) : T(T), n_t(n_t), X(X), n_x(n_x), f_a(f_a){
    h = X / n_x;
    tau = T / n_t;

    for (int n = 0; n < n_t; n++) {
        field.emplace_back(n_x);
    }

    for (int i = 0; i < n_x; i++) {
        field[0][i] = f_beg((i+0.5l)*h);
    }
}

int SolverFVM::calculate() {
    for (int n = 1; n < n_t; n++) {
        for (int i = 0; i < n_x; i++) {
            calculate_cell(n, i);
        }
    }
    return 0;
}

const std::vector<std::vector<long double>> &SolverFVM::getField() const {
    return field;
}


long double SolverGodunov::w_l(int n, int i) { // TODO: Доработать для объемов, а не узлов
    return (f_a(n * tau, i * h, field[n][i]) + f_a(n * tau, (i - 1) * h, field[n][(n_x + i - 1) % n_x])) / 2 * field[n][(n_x + i - 1) % n_x];
}

long double SolverGodunov::w_r(int n, int i) { // TODO: Доработать для объемов, а не узлов
    return (f_a(n * tau, i * h, field[n][i]) + f_a(n * tau, (i + 1) * h, field[n][(i+1) % n_x])) / 2 * field[n][i];
}

SolverGodunov::SolverGodunov(long double T, long double tau, long double X, long double h,
                             long double (*f_a)(long double, long double, long double),
                             long double (*f_beg)(long double))
        : SolverFVM(T, tau, X, h, f_a, f_beg) {

}

SolverGodunov::SolverGodunov(long double T, int n_t, long double X, int n_x,
                             long double (*f_a)(long double, long double, long double),
                             long double (*f_beg)(long double))
        : SolverFVM(T, n_t, X, n_x, f_a, f_beg) {

}
