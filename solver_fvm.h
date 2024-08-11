//
// Created by epsilon on 28.07.24.
//

#ifndef COMPMOD_SOLVER_FVM_H
#define COMPMOD_SOLVER_FVM_H


#include <functional>

class SolverFVM {
private:
    void calculate_cell(int n, int i, long double t);
protected:
    long double T1, T2, X1, X2, CFL, h;
    unsigned int n_x;
    bool isPeriodic;
    std::vector<std::vector<long double>> field;
    std::vector<long double> taus;
    
    std::function<long double(long double, long double, long double)> f_a, df_a;
    std::function<long double(long double)> f_beg, f_l, f_r;

    virtual long double w_l(int n, int i, long double t) = 0;
    virtual long double w_r(int n, int i, long double t) = 0;

public:
    SolverFVM(long double T1, long double T2, long double X1, long double X2, unsigned int n_x, long double CFL,
              long double (*f_a)(long double, long double, long double),
              long double (*df_a)(long double, long double, long double),
              long double (*f_beg)(long double));
/*
    SolverFVM(long double T1, long double T2, long double X1, long double X2, unsigned int n_x, long double CFL,
              long double(*f_a)(long double, long double, long double),
              long double(*f_l)(long double),
              long double(*f_r)(long double),
              long double(*f_beg)(long double));
*/
    int calculate();


    [[nodiscard]] const std::vector<long double> &getTaus() const;
    [[nodiscard]] const std::vector<std::vector<long double>> &getField() const;

};


class SolverGodunov : public SolverFVM {
    long double w_l(int n, int i, long double t) override;
    long double w_r(int n, int i, long double t) override;

public:
    SolverGodunov(long double T1, long double T2, long double X1, long double X2, unsigned int n_x, long double CFL,
                  long double(*f_a)(long double, long double, long double),
                  long double(*df_a)(long double, long double, long double),
                  long double(*f_beg)(long double));


};


#endif //COMPMOD_SOLVER_FVM_H
