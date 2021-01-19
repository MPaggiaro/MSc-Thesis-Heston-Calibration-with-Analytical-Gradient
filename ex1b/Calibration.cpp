//
// Created by Marco Paggiaro on 11/01/2021.
//

#include "Calibration.h"

Calibration::Calibration(std::vector<double> strikes, std::vector<double> maturities, double r, double q, double kappa,
                         double theta, double sigma, double rho, double v0, double S0) {
    options.reserve(strikes.size());
    EuropeanOption::market = Market(r,q,kappa,theta,sigma,rho,v0,S0);
    for (unsigned i = 0; i < strikes.size(); ++i)
        options.emplace_back(strikes[i], maturities[i]);
}

Calibration::Calibration(std::vector<double> strikes, std::vector<double> maturities, double *parameters, double r,
                         double S0, double q, unsigned nParameters) {
    options.reserve(strikes.size());
    for (unsigned i = 0; i < strikes.size(); ++i)
        options.emplace_back(strikes[i], maturities[i]);
    setParameters(parameters);
    EuropeanOption::market.r = r;
    EuropeanOption::market.S0 = S0;
    EuropeanOption::market.q = q;
    EuropeanOption::market.nParameters = nParameters;
}

void Calibration::setParameters(const double *parameters) {
    EuropeanOption::market.v0 = parameters[0];
    EuropeanOption::market.theta = parameters[1];
    EuropeanOption::market.rho = parameters[2];
    EuropeanOption::market.kappa = parameters[3];
    EuropeanOption::market.sigma = parameters[4];
}
