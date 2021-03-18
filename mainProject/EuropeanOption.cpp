//
// Created by Marco Paggiaro on 02/01/2021.
//

#include "EuropeanOption.h"
#include <limits>
#include <algorithm>

// definition of static members:
double EuropeanOption::q, EuropeanOption::kappa,
    EuropeanOption::theta, EuropeanOption::sigma, EuropeanOption::rho,
    EuropeanOption::v0, EuropeanOption::S0, EuropeanOption::tau_bar = 0.08219178;
unsigned EuropeanOption::nParameters;
std::vector<double> EuropeanOption::times, EuropeanOption::deltaTimes,
    EuropeanOption::phiT;
// std::vector<double> EuropeanOption::xGauss, EuropeanOption::wGauss;

EuropeanOption::EuropeanOption (const double K, const double T, const double r, std::string  optType,
                                std::string cfType,const double N):
        K(K),T(T),r(r), optType(std::move(optType)),N(N),cfType(std::move(cfType)) { }


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


