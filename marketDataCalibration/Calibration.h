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

    std::vector<double> marketPrices, calibratedPrices;
    // Parameters for the levmar algorithm:
    std::vector<double> info, opts;
    double startClock = 0.0, stopClock = 0.0;
    int maxIterations = 35;

    // parameter sets for the levmar calibration:
    std::vector<double> searchParameters;
    std::vector<double> marketParameters;
    int nParameters;

    // Dates for the Displacement vector:
    std::vector<double> optionExpiries;
    std::vector<double> deltaTimes;

    Calibration(const std::vector<double> &SPX_strikes,
                const std::vector<double> &SPX_maturities, const double *parameters, unsigned nParameters,
                double r, double S0, double q = 0.0 );

    Calibration(const std::vector<double> &SPX_strikes,  const std::vector<double> &SPX_maturities,
                const std::vector<double> &VIX_strikes, const std::vector<double> &VIX_maturities,
                double *parameters, unsigned nParameters, double r, double S0, bool displacementFlag = false,
                double q = 0.0 );
    // constructor for market prices:
    Calibration(const std::vector<double> &SPX_strikes,  const std::vector<double> &SPX_maturities,
                const std::vector<double> &SPX_mkt_prices,const std::vector<double> &SPX_rates,
                const std::vector<double> &VIX_strikes, const std::vector<double> &VIX_maturities,
                const std::vector<double> &VIX_mkt_prices,const std::vector<double> &VIX_rates,
                 const std::vector<std::string> &SPX_option_types,
                double S0, double VIX_0, bool displacementFlag = false,
                double q = 0.0 );
    // void setMarketPrices();
    void setParameters (const double *parameters);
    void setMarketParameters(const double *parameters);
    void setTermStructure(const std::vector<double> &SPX_maturities,const std::vector<double> &VIX_maturities);
    void saveCalibration (const double *parameters, const double *information);
    unsigned size() const;

    // static std::vector<double> IntegralDisplacement() ;

    std::vector<double> Prices() const;
    std::vector<double> SPX_Prices() const;
    std::vector<double> VIX_Prices() const;

    std::vector<double> Gradients() const;
    std::vector<double> SPX_Gradients() const;
    std::vector<double> VIX_Gradients() const;

    void print(const double *initialGuess) const;
};

// void computePrices(double *parameters, double *prices, int m, int n, void *data);
// void computeGradients(double *parameters, double *gradient, int m, int n, void *data);
void objectiveFunction(double *parameters, double *objectiveFunction, int m, int n, void *data);
void gradientObjective(double *parameters, double *gradient, int m, int n, void *data);
void calibrate (Calibration &calibration, const double *initialGuess, int maxIterations = 35, bool perturbation = false,
                const std::string &gradientType = "Analytical");
void calibrateMarketData (Calibration &calibration, const double *initialGuess, bool perturbation = false,
                const std::string &gradientType = "Analytical");


void testModel (Calibration &calibration, int nIterations, const std::string &gradientType = "Analytical");

// void perturbPrices(double *prices, int size);
void generateGuess(double *parameters, int size);

#endif //EX1_CALIBRATION_H
