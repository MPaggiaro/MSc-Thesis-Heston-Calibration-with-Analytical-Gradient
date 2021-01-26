//
// Created by Marco Paggiaro on 02/01/2021.
//

#ifndef EX1_EUROPEANOPTION_H
#define EX1_EUROPEANOPTION_H

#include <string>
#include <utility>
#include <vector>
#include <complex>


class EuropeanOption {
private:
    void copy(const EuropeanOption& o2); // need to check if I can do this with inheritance, see book.

public:
// Public member data for convenience only
//    static Market market;
    static double r, q, kappa, theta, sigma, rho, v0, S0;
    static unsigned nParameters;
    double K = 0; // Strike price
    double T = 0; // Expiry date
    std::string optType; // Option name (call, put)

    double N = 200;
    std::string cfType;
//    static std::vector<double> xGauss, wGauss;

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
    virtual double Price() const = 0; // some prices are negative, have to fix this!
    // virtual double ImprovedPrice() const = 0;
    virtual std::vector<double> Jacobian() const = 0;

private:
    virtual std::complex<double> CharFunc(std::complex<double> u) const = 0;
    virtual std::vector<std::complex<double>> JacobianCF(std::complex<double> u) const = 0;
};

#endif //EX1_EUROPEANOPTION_H
