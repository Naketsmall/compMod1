//
// Created by epsilon on 28.07.24.
//

#include "verifier.h"
#include "cmath"

Verifier::Verifier(std::vector<std::vector<long double>> &calculated,
                   long double T, long double tau, long double X, long double h, long double a,
                   long double(*f_beg)(long double))
: calculated(calculated), tau(tau), h(h) {
    accurate = calculated;
    for (int n = 1; n < calculated.size(); n++) {
        for (int i = 0; i < calculated[0].size(); i++) {
            accurate[n][i] = f_beg((i*h - n*tau*a) -  (int)((i*h - n*tau*a)/X)*X);
        }
    }

    diff = calculated;
    for (int n = 0; n < calculated.size(); n++) {
        for (int i = 0; i < calculated[0].size(); i++) {
            diff[n][i] = calculated[n][i] - accurate[n][i];
        }
    }

}

const std::vector<std::vector<long double>> &Verifier::getAccurate() const {
    return accurate;
}

const std::vector<std::vector<long double>> &Verifier::getDiff() const {
    return diff;
}

long double Verifier::get_max() {
    long double max = 0;
    for (int n = 0; n < calculated.size(); n++) {
        for (int i = 0; i < calculated[0].size(); i++) {
            max = std::abs(diff[n][i]) > max ? std::abs(diff[n][i]) : max;
        }
    }
    return max;
}

long double Verifier::get_L1() {
    long double sum = 0;
    for (int n = 0; n < calculated.size(); n++) {
        for (int i = 0; i < calculated[0].size(); i++) {
            sum += diff[n][i];
        }
    }
    return sum*h;
}
