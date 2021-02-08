//
// Created by Marco Paggiaro on 11/01/2021.
//

#ifndef EX1_CALIBRATION_H
#define EX1_CALIBRATION_H
#include <utility>
#include <vector>
#include <memory>
#include "EuropeanOption.h"
#include "SPXOption.h"
#include "VIXOption.h"

class Calibration {
public:
    std::vector<SPXOption> SPX_options;
    std::vector<VIXOption> VIX_options;

    Calibration(const std::vector<double> &SPX_strikes,
                const std::vector<double> &SPX_maturities, double *parameters,
                double r, double S0, double q = 0.0, unsigned nParameters = 5);

    Calibration(const std::vector<double> &SPX_strikes,  const std::vector<double> &SPX_maturities,
                const std::vector<double> &VIX_strikes, const std::vector<double> &VIX_maturities,
                double *parameters, double r, double S0, bool displacementFlag = false,
                double q = 0.0, unsigned nParameters = 5);

    static void setParameters (const double *parameters);
    unsigned size() const;

    std::vector<double> Prices() const;
    std::vector<double> SPX_Prices() const;
    std::vector<double> VIX_Prices() const;

    std::vector<double> Gradients() const;
    std::vector<double> SPX_Gradients() const;
    std::vector<double> VIX_Gradients() const;
};

void computePrices(double *parameters, double *prices, int m, int n, void *data);
void computeGradients(double *parameters, double *gradient, int m, int n, void *data);
void calibrate (const Calibration &calibration, const double *initialGuess,
                const std::string &gradientType = "Analytical");


#endif //EX1_CALIBRATION_H
