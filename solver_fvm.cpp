//
// Created by epsilon on 28.07.24.
//

#include <cmath>
#include "solver_fvm.h"
#include "stdio.h"

/*
SolverFVM::SolverFVM(long double T1, long double T2, long double X1, long double X2, unsigned int n_x, long double CFL,
                     long double (*f_a)(long double, long double, long double), long double (*f_l)(long double),
                     long double (*f_r)(long double), long double (*f_beg)(long double))
                     : T1(T1), T2(T2), X1(X1), X2(X2), n_x(n_x), CFL(CFL), f_a(f_a), f_beg(f_beg), f_l(f_l), f_r(f_r) {
    h = (X2 - X1) / n_x;
    taus = {0};
    isPeriodic = false;
}
 */

SolverFVM::SolverFVM(long double T1, long double T2, long double X1, long double X2, unsigned int n_x, long double CFL,
                     long double (*f_a)(long double, long double, long double),
                     long double (*df_a)(long double, long double, long double),
                     long double (*f_beg)(long double))
                     : T1(T1), T2(T2), X1(X1), X2(X2), n_x(n_x), CFL(CFL), f_a(f_a), df_a(df_a), f_beg(f_beg) {
    h = (X2 - X1) / n_x;
    taus = {0};
    isPeriodic = true;
}

void SolverFVM::calculate_cell(int n, int i, long double t) { // works for layers 1, 2... not for 0
    field[n][i] = field[n - 1][i] - taus[n] / h * (w_r(n - 1, i, t) - w_l(n - 1, i, t));
}



const std::vector<std::vector<long double>> &SolverFVM::getField() const {
    return field;
}

int SolverFVM::calculate() { // Посчитано для периодичных ГУ
    long double min_tau = 1000000000000l;
    long double t = T1;
    bool last_iteration = false;
    field.emplace_back(n_x);
    for (int i = 0; i < n_x; i++) {
        field[0][i] = f_beg(X1+(i + 0.5l) * h);
        min_tau = std::fmin(min_tau, std::abs(CFL * h / df_a(T1, X1+(i + 0.5l) * h, field[0][i])));
    }
    t += min_tau;
    taus.push_back(min_tau);
    min_tau = 1000000000000l;

    for (int n = 1; n > 0; n++) {
        field.emplace_back(n_x);
        for (int i = 0; i < n_x; i++) {
            calculate_cell(n, i, t);
            min_tau = std::fmin(min_tau, std::abs(CFL * h / df_a(t, X1+(i + 0.5l) * h, field[n][i])));
        }

        if (last_iteration)
            break;
        if (t + min_tau > T2) {
            last_iteration = true;
            min_tau = T2 - t;
        }
        t += min_tau;
        //printf("t: %Lf, min_tau: %Lf, field[%d][0]: %Lf, last_iteration %b\n", t, min_tau, n, f_a(t, X1+(0 + 0.5l) * h, field[n][50]), last_iteration);
        taus.push_back(min_tau);
        min_tau = 1000000000000l;
    }
    return 0;
}

const std::vector<long double> &SolverFVM::getTaus() const {
    return taus;
}


long double SolverGodunov::w_l(int n, int i, long double t) { // TODO: Доработать для объемов, а не узлов
    //return f_a(t, X1+i*h, (field[n][(n_x + i - 1) % n_x]+field[n][i])/2);
    if (df_a(taus[n], i * h, field[n][i]) + df_a(taus[n], (i - 1) * h, field[n][(n_x + i - 1) % n_x]) >= 0)
        return (df_a(taus[n], i * h, field[n][i]) + df_a(taus[n], (i - 1) * h, field[n][(n_x + i - 1) % n_x])) / 2 *
            field[n][(n_x + i - 1) % n_x];
    else
        return (df_a(taus[n], i * h, field[n][i]) + df_a(taus[n], (i - 1) * h, field[n][(n_x + i - 1) % n_x])) / 2 *
               field[n][i];
}

long double SolverGodunov::w_r(int n, int i, long double t) { // TODO: Доработать для объемов, а не узлов
    //return f_a(t, X1+(i+1)*h, (field[n][i]+field[n][(i+1) % n_x])/2);
    if (df_a(taus[n], i * h, field[n][i]) + df_a(taus[n], (i + 1) * h, field[n][(i+1) % n_x]) >= 0)
        return (df_a(taus[n], i * h, field[n][i]) + df_a(taus[n], (i + 1) * h, field[n][(i+1) % n_x])) / 2 *
            field[n][i];
    else
        return (df_a(taus[n], i * h, field[n][i]) + df_a(taus[n], (i + 1) * h, field[n][(i+1) % n_x])) / 2 *
               field[n][(n_x + i + 1) % n_x];
}
/*
long double SolverGodunov::w_l(int n, int i, long double t) {
    if (df_a(t, X1 + i * h, (field[n][i] + field[n][(n_x+i-1) % n_x])/2) >= 0) {
        printf("w_l: bigger\n");
        return f_a(t, X1 + i * h, field[n][(n_x + i - 1) / n_x]);
    }
    else
        return f_a(t, X1 + i * h, field[n][i]);
}

long double SolverGodunov::w_r(int n, int i, long double t) {
    if (df_a(t, X1 + (i+1) * h, (field[n][i] + field[n][(i+1) % n_x])/2) >= 0) {
        printf("w_r: bigger\n");
        return f_a(t, X1 + (i + 1) * h, field[n][i]);
    }
    else
        return f_a(t, X1 + (i+1) * h, field[n][(i+1) % n_x]);
}*/

SolverGodunov::SolverGodunov(long double T1, long double T2, long double X1, long double X2, unsigned int n_x,
                             long double CFL,
                             long double (*f_a)(long double, long double, long double),
                             long double (*df_a)(long double, long double, long double),
                             long double (*f_beg)(long double)) : SolverFVM(T1, T2, X1, X2, n_x, CFL, f_a, df_a, f_beg) {

}


