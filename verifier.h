//
// Created by epsilon on 28.07.24.
//

#ifndef COMPMOD_VERIFIER_H
#define COMPMOD_VERIFIER_H

#include <vector>

/// должен сравнивать точное и посчитанное решения, находить нормы max, L1 от разницы решений
///
class Verifier {
private:
    long double T, tau, X, h;
    std::vector<std::vector<long double>> accurate;
    std::vector<std::vector<long double>> calculated;
    std::vector<std::vector<long double>> diff;

public:
    Verifier(std::vector<std::vector<long double>> &calculated,
             long double T, long double tau, long double X, long double h, long double a,
             long double(*f_beg)(long double));

    [[nodiscard]] const std::vector<std::vector<long double>> &getAccurate() const;
    [[nodiscard]] const std::vector<std::vector<long double>> &getDiff() const;

    long double get_max();
    long double get_L1();



};


#endif //COMPMOD_VERIFIER_H
