//
// Created by Marco Paggiaro on 19/01/2021.
//

#include "VIXOption.h"
#include <complex>

double VIXOption::Price() const {
    return 0;
}

std::vector<double> VIXOption::Jacobian() const {
    return std::vector<double>();
}

std::complex<double> VIXOption::CharFunc(std::complex<double> u) const {
    // Renaming the market variables for better readability.
    auto &r = market.r, &q = market.q, &S0 = market.S0, &kappa = market.kappa,
            &sigma = market.sigma, &rho = market.rho, &theta = market.theta, &v0 = market.v0;

    double sigma2 = pow(sigma,2);
    double ekt = exp(-kappa*T);
    return exp(-2*kappa*theta/sigma2 * log(1.0 - I*u*sigma2/(2*kappa)*(1 - ekt))
            + v0 * (I*u*ekt)/(1.0 - I*u*sigma2/(2*kappa)*(1 - ekt)));
}

std::vector<std::complex<double>> VIXOption::JacCharFunc(std::complex<double> u) const {
    return std::vector<std::complex<double>>();
}
