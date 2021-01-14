//
// Created by Marco Paggiaro on 11/01/2021.
//

#include "Calibration.h"

//std::vector<double> Calibration::GetPrices() const {
//    std::vector<double> prices(options.size());
//    for (unsigned i = 0; i < options.size(); ++i)
//    {
//        prices[i] = options[i].Price();
//    }
//    return prices;
//}

//std::vector<double> Calibration::GetGradients() const {
//
//    // For the moment, we will return a vector of size M * N, where:
//    // M = number of market parameters, N = number of observations.
//    const unsigned M = EuropeanOption::market.nParameters;
//    std::vector<double> gradients(options.size() * M);
//    for (unsigned i = 0; i < options.size(); ++i)
//    {
//        auto currentGradient = options[i].Jacobian();
//        for (unsigned j = 0; j < M; ++j)
//        {
//            gradients[M * i + j] = currentGradient[j];
//        }
//    }
//    return gradients;
//}

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

//
//void Calibration::Calibrate() const {
//    double marketPrices[options.size()];
//
//    // Getting the parameters vector.
//    double parameters[EuropeanOption::market.nParameters];
//    parameters[0] = EuropeanOption::market.v0;
//    parameters[1] = EuropeanOption::market.theta;
//    parameters[2] = EuropeanOption::market.rho;
//    parameters[3] = EuropeanOption::market.kappa;
//    parameters[4] = EuropeanOption::market.sigma;
//
//    fHes(parameters,marketPrices,EuropeanOption::market.nParameters,options.size(),
//         nullptr);
//
//    // algorithm parameters
//    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
//    opts[0]=LM_INIT_MU;
//    // stopping thresholds for
//    opts[1]=1E-10;       // ||J^T e||_inf
//    opts[2]=1E-10;       // ||Dp||_2
//    opts[3]=1E-10;       // ||e||_2
//    opts[4]= LM_DIFF_DELTA; // finite difference if used
//
//    double initialGuess[EuropeanOption::market.nParameters];
//    initialGuess[0] = 0.2;
//    initialGuess[1] = 0.2;
//    initialGuess[2] = -0.6;
//    initialGuess[3] = 1.2;
//    initialGuess[4] = 0.3;
//
//    dlevmar_der(fHes,JacHes,initialGuess,marketPrices,EuropeanOption::market.nParameters,
//                options.size(),100,opts,info, NULL, NULL, NULL);
//
//}
