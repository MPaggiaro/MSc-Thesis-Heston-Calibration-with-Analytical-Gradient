//
// Created by Marco Paggiaro on 11/01/2021.
//

#include "Calibration.h"
#include "SPXOption.h"
//#include <boost/math/quadrature/gauss.hpp>

Calibration::Calibration(std::vector<double> strikes, std::vector<double> maturities, double *parameters, double r,
                         double S0, double q, unsigned nParameters) {
    options.reserve(strikes.size());
    for (unsigned i = 0; i < strikes.size(); ++i)
        options.emplace_back(
                std::make_shared<SPXOption>(strikes[i], maturities[i]));
    setParameters(parameters);
    EuropeanOption::r = r;
    EuropeanOption::S0 = S0;
    EuropeanOption::q = q;
    EuropeanOption::nParameters = nParameters;

//    EuropeanOption::xGauss = boost::math::quadrature::gauss<double, 64>::abscissa();
//    EuropeanOption::wGauss = boost::math::quadrature::gauss<double, 64>::weights();

}

void Calibration::setParameters(const double *parameters) {
    EuropeanOption::v0 = parameters[0];
    EuropeanOption::theta = parameters[1];
    EuropeanOption::rho = parameters[2];
    EuropeanOption::kappa = parameters[3];
    EuropeanOption::sigma = parameters[4];
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
