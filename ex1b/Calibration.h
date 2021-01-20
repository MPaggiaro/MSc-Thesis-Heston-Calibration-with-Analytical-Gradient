//
// Created by Marco Paggiaro on 11/01/2021.
//

#ifndef EX1_CALIBRATION_H
#define EX1_CALIBRATION_H
#include <utility>
#include <vector>
#include <memory>
#include "Market.h"
#include "EuropeanOption.h"

class Calibration {
public:
    std::vector<std::shared_ptr<EuropeanOption>> options;

    Calibration(std::vector<double> strikes, std::vector<double> maturities, double *parameters,
                double r, double S0, double q = 0.0, unsigned nParameters = 5);

    static void setParameters (const double *parameters);

};

void computePrices(double *p, double *x, int m, int n, void *data);
void computeGradients(double *p, double *jac, int m, int n, void *data);


#endif //EX1_CALIBRATION_H
