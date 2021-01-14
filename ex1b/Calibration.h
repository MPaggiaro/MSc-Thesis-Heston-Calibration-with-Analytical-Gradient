//
// Created by Marco Paggiaro on 11/01/2021.
//

#ifndef EX1_CALIBRATION_H
#define EX1_CALIBRATION_H
#include <utility>
#include <vector>
#include "Market.h"
#include "EuropeanOption.h"

class Calibration {
public:
    std::vector<EuropeanOption> options;

    Calibration(std::vector<double> strikes, std::vector<double> maturities,
                double r, double q, double kappa, double theta, double sigma, double rho,
                double v0, double S0);
    Calibration(std::vector<double> strikes, std::vector<double> maturities, double *parameters,
                double r, double S0, double q = 0.0, unsigned nParameters = 5);

    // std::vector<double> GetPrices () const;
    // std::vector<double> GetGradients () const;
    static void setParameters (const double *parameters);

private:

};


#endif //EX1_CALIBRATION_H
