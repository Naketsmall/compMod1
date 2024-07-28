#include <iostream>
#include <cmath>
#include "solver_fvm.h"

long double f_a(long double t, long double x) {
    return 1;
}

long double f_beg(long double x) {
    return sin(2*3.1415*x);
}

int main() {
    SolverGodunov s1 = SolverGodunov(1, 0.1, 1, 0.1, &f_a, &f_beg);
    s1.calculate();
    printf("%d\n", 8-1 % 8);

    std::vector<std::vector<long double>> field = s1.getField();
    for (int n = 0; n < field.size(); n++) {
        for (int i = 0; i < field[0].size(); i++) {
            printf("%Lf ", field[n][i]);
        }
        printf("\n");
    }


}
