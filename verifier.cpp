//
// Created by epsilon on 28.07.24.
//

#include <cstdio>
#include <sstream>
#include <iomanip>
#include "verifier.h"
#include "cmath"

long double ld_fmod(long double a, long double b) {
    return a - static_cast<long long>(a / b) * b;
}

long long ld_fdiv(long double a, long double b) {
    return static_cast<long long>((a) / b);
}

Verifier::Verifier(std::vector<std::vector<long double>> &calculated,
                   long double T, long double tau, long double X, long double h, long double a,
                   long double(*f_beg)(long double))
: calculated(calculated), T(T), tau(tau), X(X), h(h) {
    accurate = calculated;
    for (int n = 1; n < calculated.size(); n++) {
        for (int i = 0; i < calculated[0].size(); i++) {
            accurate[n][i] = f_beg(((i+0.5l)*h - n*tau*a) +
                    ((i+0.5l)*h-n*tau*a < 0 ? X * (ld_fdiv(-((i+0.5l)*h-n*tau*a), X) + 1): 0));
              //accurate[n][i] = f_beg(ld_fmod((i*h - n*tau*a), X));
        }
    }
    //L1: 0.00008705738838428928
    //Max: 0.00018530717932098052

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
            sum += std::abs(diff[n][i]);
        }
    }
    return sum / std::ceil(X / h) / std::ceil(T / tau);
}

std::string vector2d_to_string(const std::vector<std::vector<long double>> &field) {
    std::stringstream line;

    for (int n = 0; n < field.size(); n++) {
        for (int i = 0; i < field[0].size(); i++) {
            if (i != field[0].size()-1)
                line << std::fixed << std::setprecision(30) << field[n][i] << ",";// field[n][i] + ",";
            else
                line << std::fixed << std::to_string(field[n][i]);
        }
        line << "\n";
    }
    return std::move(line).str();
}