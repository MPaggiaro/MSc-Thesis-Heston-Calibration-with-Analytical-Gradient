//
// Created by Marco Paggiaro on 02/01/2021.
//

#include "EuropeanOption.h"

Market EuropeanOption::market;

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


