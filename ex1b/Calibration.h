//
// Created by Marco Paggiaro on 11/01/2021.
//

#ifndef EX1_CALIBRATION_H
#define EX1_CALIBRATION_H
#include <utility>
#include <vector>
#include "Market.h"

class Calibration {
public:
    std::vector<double> strikes;
    std::vector<double> maturities;
    const Market market;

    Calibration(std::vector<double> strikes, std::vector<double> maturities,
                double r, double q, double kappa, double theta, double sigma, double rho,
                double v0, double S0): strikes(std::move(strikes)), maturities(std::move(maturities)),
                market(r,q,kappa,theta,sigma,rho,v0,S0) {}
};


#endif //EX1_CALIBRATION_H
