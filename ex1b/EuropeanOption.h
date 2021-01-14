//
// Created by Marco Paggiaro on 02/01/2021.
//

#ifndef EX1_EUROPEANOPTION_H
#define EX1_EUROPEANOPTION_H

#include <string>
#include <utility>
#include <complex>
#include "Market.h"

class EuropeanOption {
private:
    void copy(const EuropeanOption& o2);

public:
// Public member data for convenience only
    // static double r, q, kappa, theta, sigma, rho, v0, S0;
    static Market market;
    double K = 0; // Strike price
    double T = 0; // Expiry date
    std::string optType; // Option name (call, put)

    // integration and CF parameters:
    static constexpr unsigned int M = 64;
    double N = 200;
    std::string cfType;

public:
// Constructors
    EuropeanOption(const double K, const double T, std::string  optType = "C",
                   std::string cfType = "Cui",const double N = 200):
                   K(K),T(T),optType(std::move(optType)),N(N),cfType(std::move(cfType)){}

// Destructor
    virtual ~EuropeanOption();
// Assignment operator
    EuropeanOption& operator = (const EuropeanOption& option2);
// Functions that calculate option price and (some) sensitivities
    double Price() const; // some prices are negative, have to fix this!
    std::vector<double> Jacobian() const;

private:
    std::complex<double> CharFunc(std::complex<double> u) const;
    std::vector<std::complex<double>> JacCharFunc(std::complex<double> u) const;
};

#endif //EX1_EUROPEANOPTION_H
