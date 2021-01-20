//
// Created by Marco Paggiaro on 11/01/2021.
//

#include "Calibration.h"
#include "SPXOption.h"

Calibration::Calibration(std::vector<double> strikes, std::vector<double> maturities, double *parameters, double r,
                         double S0, double q, unsigned nParameters) {
    options.reserve(strikes.size());
    for (unsigned i = 0; i < strikes.size(); ++i)
        options.emplace_back(
                std::make_shared<SPXOption>(strikes[i], maturities[i]));
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

void computePrices(double *p, double *x, int m, int n, void *data)
{
    auto * calibration = static_cast<Calibration *> (data);
    Calibration::setParameters(p);
    for (unsigned i = 0; i < calibration->options.size(); ++i)
    {
        x[i] = calibration->options[i]->Price();
    }
}

void computeGradients(double *p, double *jac, int m, int n, void *data)
{
    auto * calibration = static_cast<Calibration *> (data);
    Calibration::setParameters(p);
    for (unsigned i = 0; i < n; ++i)
    {
        auto currentGradient = calibration->options[i]->Jacobian();
        for (unsigned j = 0; j < m; ++j)
        {
            jac[m * i + j] = currentGradient[j];
        }
    }
}
