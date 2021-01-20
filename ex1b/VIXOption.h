//
// Created by Marco Paggiaro on 19/01/2021.
//

#ifndef EX1_VIXOPTION_H
#define EX1_VIXOPTION_H
#include <utility>

#include "EuropeanOption.h"

class VIXOption : public EuropeanOption {
public:
    // Constructor:
    VIXOption(const double K, const double T, std::string  optType = "C",
              std::string cfType = "Cui",const double N = 200):
            EuropeanOption(K,T,std::move(optType),std::move(cfType),N){ }

    double Price() const override;
    std::vector<double> Jacobian() const override;

private:
    std::complex<double> CharFunc(std::complex<double> u) const override;
    std::vector<std::complex<double>> JacCharFunc(std::complex<double> u) const override;
};


#endif //EX1_VIXOPTION_H
