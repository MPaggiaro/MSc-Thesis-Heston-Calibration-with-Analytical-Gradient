//
// Created by Marco Paggiaro on 11/01/2021.
//

#ifndef EX1_CALIBRATION_H
#define EX1_CALIBRATION_H
#include <utility>
#include <vector>
#include <memory>
#include <levmar/levmar.h>
#include "EuropeanOption.h"
#include "SPXOption.h"
#include "VIXOption.h"

class Calibration {
public:
    std::vector<SPXOption> SPX_options;
    std::vector<VIXOption> VIX_options;

    std::vector<double> marketPrices;
    // Parameters for the levmar algorithm:
    std::vector<double> info, opts;
    double startClock = 0.0, stopClock = 0.0;

    // parameter sets for the levmar calibration:
    std::vector<double> searchParameters;

    Calibration(const std::vector<double> &SPX_strikes,
                const std::vector<double> &SPX_maturities, double *parameters,
                double r, double S0, double q = 0.0, unsigned nParameters = 5);

    Calibration(const std::vector<double> &SPX_strikes,  const std::vector<double> &SPX_maturities,
                const std::vector<double> &VIX_strikes, const std::vector<double> &VIX_maturities,
                double *parameters, double r, double S0, bool displacementFlag = false,
                double q = 0.0, unsigned nParameters = 5);

    void setMarketPrices();
    static void setParameters (const double *parameters);
    void saveCalibration (const double *parameters, const double *information);
    unsigned size() const;

    static std::vector<double> IntegralDisplacement() ;

    std::vector<double> Prices() const;
    std::vector<double> SPX_Prices() const;
    std::vector<double> VIX_Prices() const;

    std::vector<double> Gradients() const;
    std::vector<double> SPX_Gradients() const;
    std::vector<double> VIX_Gradients() const;

    void print(const double *initialGuess) const;
};

void computePrices(double *parameters, double *prices, int m, int n, void *data);
void computeGradients(double *parameters, double *gradient, int m, int n, void *data);
void objectiveFunction(double *parameters, double *objectiveFunction, int m, int n, void *data);
void gradientObjective(double *parameters, double *gradient, int m, int n, void *data);
void calibrate (Calibration &calibration, const double *initialGuess, bool perturbation = false,
                const std::string &gradientType = "Analytical");
void testModel (Calibration &calibration, int nIterations, const std::string &gradientType = "Analytical");

void perturbPrices(double *prices, int size);
void generateGuess(double *parameters, int size);

#endif //EX1_CALIBRATION_H
