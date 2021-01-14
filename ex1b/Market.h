//
// Created by Marco Paggiaro on 11/01/2021.
//

#ifndef EX1_MARKET_H
#define EX1_MARKET_H


class Market {
public:
    // Market parameters:
    double r, q, kappa, theta, sigma, rho, v0, S0;

    // number of unknown parameters:
    unsigned nParameters;

    Market() = default;
    Market(double r, double q, double kappa, double theta, double sigma, double rho, double v0,
           double S0, unsigned nParameters = 5): r(r), q(q), kappa(kappa), theta(theta), sigma(sigma), rho(rho), v0(v0),
           S0(S0), nParameters(nParameters) {}
};


#endif //EX1_MARKET_H
