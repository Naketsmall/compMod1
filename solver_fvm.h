//
// Created by epsilon on 28.07.24.
//

#ifndef COMPMOD_SOLVER_FVM_H
#define COMPMOD_SOLVER_FVM_H


#include <functional>

class SolverFVM {
private:
    void calculate_cell(int n, int i);
protected:
    long double T, tau, X, h;
    unsigned int n_t, n_x;
    std::vector<std::vector<long double>> field;

    std::function<long double(long double, long double, long double)> f_a;

    virtual long double w_l(int n, int i) = 0;
    virtual long double w_r(int n, int i) = 0;

public:
    SolverFVM(long double T, long double tau, long double X, long double h,
              long double(*f_a)(long double, long double, long double),
              long double(*f_beg)(long double));

    SolverFVM(long double T, int n_t, long double X, int n_x,
              long double(*f_a)(long double, long double, long double),
              long double(*f_beg)(long double));

    int calculate();

    long double get_courant(long double t, long double x, long double u);

    [[nodiscard]] const std::vector<std::vector<long double>> &getField() const;

};

class SolverGodunov : public SolverFVM {
    long double w_l(int n, int i) override;
    long double w_r(int n, int i) override;

public:
    SolverGodunov(long double T, long double tau, long double X, long double h,
                  long double(*f_a)(long double, long double, long double),
                  long double(*f_beg)(long double));

    SolverGodunov(long double T, int n_t, long double X, int n_x,
                  long double(*f_a)(long double, long double, long double),
                  long double(*f_beg)(long double));

};


#endif //COMPMOD_SOLVER_FVM_H
