//
// Created by Marco Paggiaro on 02/01/2021.
//

#include "EuropeanOption.h"

//Market EuropeanOption::market;

double EuropeanOption::r, EuropeanOption::q, EuropeanOption::kappa,
    EuropeanOption::theta, EuropeanOption::sigma, EuropeanOption::rho,
    EuropeanOption::v0, EuropeanOption::S0;
unsigned EuropeanOption::nParameters;
// std::vector<double> EuropeanOption::xGauss, EuropeanOption::wGauss;

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


