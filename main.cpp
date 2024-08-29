#include <iostream>
#include <fstream>
#include <iomanip>
#include "solver_fvm.h"
#include "verifier.h"
#include "input.h"


int main() {

    SolverGodunov s1 = SolverGodunov(T1, T2, X1, X2, n_x, CFL, &f_a, &df_a, &f_beg_disc);
    s1.calculate();
    std::vector<std::vector<long double>> field = s1.getField();

    //std::ofstream out_acc(outExactFilename);
    //out_acc << vector2d_to_string(v1.getAccurate());
    //out_acc.close();

    std::ofstream out_calc(outCalcFilename);
    out_calc << vector2d_to_string(field);
    out_calc.close();

    std::vector<std::vector<long double>> taus2;
    taus2.push_back(s1.getTaus());
    std::ofstream out_taus("taus.csv");
    out_taus << vector2d_to_string(taus2);
    out_taus.close();
    //std::ofstream out_diff(outDiffFilename);
    //out_diff << vector2d_to_string(v1.getDiff());
    //out_diff.close();

    //printf("Max: %.30Lf\n", v1.get_max());
    //printf("L1: %.30Lf\n", v1.get_L1());
    printf("shape calc t:%zu, x:%zu\n", field.size(), field[0].size());
    //printf("shape acc: %zu, %zu\n", v1.getAccurate().size(), v1.getAccurate()[0].size());


}
