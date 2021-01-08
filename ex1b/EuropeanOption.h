//
// Created by Marco Paggiaro on 02/01/2021.
//

#ifndef EX1_EUROPEANOPTION_H
#define EX1_EUROPEANOPTION_H

#include <string>
#include <utility>
#include <complex>

class EuropeanOption {
private:
    void copy(const EuropeanOption& o2);

public:
// Public member data for convenience only
    double r = 0; // Interest rate
    double kappa = 0;
    double theta = 0;
    double sigma = 0;
    double rho = 0;
    double v0 = 0;
    double K = 0; // Strike price
    double T = 0; // Expiry date
    double S0 = 0; // Current underlying price
    double q = 0;
    std::string optType; // Option name (call, put)

    // integration and CF parameters:
    static constexpr unsigned int M = 64;
    double N = 200;
    std::string cfType;

public:
// Constructors
    EuropeanOption(const double r, const double kappa, const double theta,
                   const double sigma, const double rho, const double v0,
                   const double K, const double T, const double S0,
                   const double q, std::string  optType,
                   std::string cfType = "Cui",const double N = 200):
                   r(r), kappa(kappa), theta(theta),sigma(sigma),rho(rho), v0(v0),K(K),T(T),
                   S0(S0),q(q),optType(std::move(optType)),N(N),cfType(std::move(cfType)){}

// Destructor
    virtual ~EuropeanOption();
// Assignment operator
    EuropeanOption& operator = (const EuropeanOption& option2);
// Functions that calculate option price and (some) sensitivities
    double Price() const;
    std::vector<double> ComputeJacobian() const;
private:
    std::complex<double> CharFunc(std::complex<double> u) const;
    std::vector<std::complex<double>> JacCharFunc(std::complex<double> u) const;
};


#endif //EX1_EUROPEANOPTION_H
