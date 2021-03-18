//
// Created by Marco Paggiaro on 02/01/2021.
//

#include "EuropeanOption.h"
#include <limits>
#include <algorithm>

// definition of static members:
double EuropeanOption::tau_bar = 0.08219178;
// std::vector<double> EuropeanOption::xGauss, EuropeanOption::wGauss;

EuropeanOption::EuropeanOption (const double K, const double T, const double r, double kappa, double theta,
                                double sigma, double rho, double v0, double S0, std::string  optType,
                                unsigned nParameters,
                                double q,
                                std::string cfType,const double N):
        K(K),T(T),r(r),
        kappa(kappa), theta(theta), sigma(sigma), rho(rho), v0(v0), S0(S0), q(q), nParameters(nParameters),
        optType(std::move(optType)),N(N),cfType(std::move(cfType)) {}

EuropeanOption::EuropeanOption(double K, double T, double r, double S0, std::string optType, unsigned int nParameters,
                               double q, std::string cfType, double N):
        K(K), T(T), r(r), S0(S0), optType(std::move(optType)),N(N),cfType(std::move(cfType)),
        q(q), nParameters(nParameters){}



void EuropeanOption::copy(const EuropeanOption& o2)
{
    K = o2.K;
    T = o2.T;
}

EuropeanOption::~EuropeanOption() = default;

EuropeanOption& EuropeanOption::operator = (const EuropeanOption& opt2)
{ // Assignment operator (deep copy)
    if (this == &opt2) return *this;
    copy (opt2);
    return *this;
}

void EuropeanOption::SetIndexT()
{
    // find i such that times[i] = T:
    auto it = std::find_if(times.begin(),times.end(),
                           [&](double time)
                           { return fabs(T - time) < std::numeric_limits<double>::epsilon(); });
    indexT = it - times.begin();
}

void EuropeanOption::setParameters(const double *parameters)
{
    v0 = parameters[0];
    theta = parameters[1];
    rho = parameters[2];
    kappa = parameters[3];
    sigma = parameters[4];

    if (nParameters > 5)
        // If we have to specify also the displacement, then we save it from here:
        for (int i = 5; i < nParameters; ++i)
            phiT[i - 5] = parameters[i];
}

void EuropeanOption::setDisplacement(const double *parameters, const std::vector<double> &finalTimes,
                                     const std::vector<double> &deltaT)
{
    times = finalTimes;
    deltaTimes = deltaT;
    phiT.resize(deltaTimes.size());
    if (parameters)
        for (int i = 5; i < nParameters; ++i)
        {
            phiT[i - 5] = parameters[i];
            // Watch out: we haven't placed (until now) a check on the size of parameters.
            // You have to place a correct number of parameters.
        }
}



