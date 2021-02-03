//
// Created by Marco Paggiaro on 11/01/2021.
//

#include <iostream>
#include <iomanip>
#include <levmar/levmar.h>
#include <algorithm>
#include "Calibration.h"
//#include <boost/math/quadrature/gauss.hpp>

Calibration::Calibration(std::vector<double> strikes, std::vector<double> maturities, double *parameters,
                         double r, double S0, double q, unsigned nParameters) {
    SPX_options.reserve(strikes.size());
    for (unsigned i = 0; i < strikes.size(); ++i)
        SPX_options.emplace_back(strikes[i], maturities[i]);
    setParameters(parameters);
    EuropeanOption::r = r;
    EuropeanOption::S0 = S0;
    EuropeanOption::q = q;
    EuropeanOption::nParameters = nParameters;
}

Calibration::Calibration(std::vector<double> SPX_strikes, std::vector<double> SPX_maturities,
                         std::vector<double> VIX_strikes, std::vector<double> VIX_maturities, double *parameters,
                         double r, double S0, double q, unsigned int nParameters):
        Calibration(std::move(SPX_strikes), std::move(SPX_maturities), parameters, r, S0, q, nParameters) {
    VIX_options.reserve(VIX_strikes.size());
    for (unsigned i = 0; i < VIX_strikes.size(); ++i)
        VIX_options.emplace_back(VIX_strikes[i], VIX_maturities[i]);
}

Calibration::Calibration(std::vector<double> SPX_strikes, const std::vector<double> &SPX_maturities,
                         std::vector<double> VIX_strikes, const std::vector<double> &VIX_maturities,
                         std::vector<double> displacement_vector, double *parameters, double r, double S0, double q,
                         unsigned int nParameters):
                         Calibration(std::move(SPX_strikes), SPX_maturities,
                                     std::move(VIX_strikes), VIX_maturities,
                                     parameters, r, S0, q, nParameters){
    auto SPX_t(SPX_maturities), VIX_t(VIX_maturities);
    SPX_t.push_back(0.0);
    std::sort( SPX_t.begin(), SPX_t.end() );
    SPX_t.erase( std::unique(SPX_t.begin(), SPX_t.end()), SPX_t.end());
    std::sort( VIX_t.begin(), VIX_t.end());
    VIX_t.erase( std::unique(VIX_t.begin(), VIX_t.end()), VIX_t.end());
    std::vector<double> partialT(SPX_t.size() + VIX_t.size());
    auto it = std::set_union(SPX_t.begin(), SPX_t.end(),
                             VIX_t.begin(), VIX_t.end(), partialT.begin());
    partialT.resize(it - partialT.begin());

    auto VIX_taubar(VIX_t);
    for (double & i : VIX_taubar) {
        i += EuropeanOption::tau_bar;
    }

    std::vector<double> finalT(partialT.size() + VIX_taubar.size());
    // unite also B + tau_bar:
    auto it2 = std::set_union(partialT.begin(), partialT.end(), VIX_taubar.begin(),
                              VIX_taubar.end(), finalT.begin());
    finalT.resize(it2 - finalT.begin());

// We need to implement a function that compares two very similar doubles:
    std::vector<int> areSame;
    for (int i = 0; i < finalT.size()-1; ++i) {
        if (AreSame(finalT[i], finalT[i+1]))
            areSame.push_back(i);
    }
    for (int i : areSame) {
        finalT.erase(finalT.begin() + i + 1);
    }

    // compute differences:
    std::vector<double> deltaT(finalT.size() - 1);
    for (int i = 0; i < deltaT.size(); ++i) {
        deltaT[i] = finalT[i+1] - finalT[i];
    }
    EuropeanOption::times = finalT;
    EuropeanOption::deltaTimes = deltaT;
    EuropeanOption::phiT = std::move(displacement_vector);
}

void Calibration::setParameters(const double *parameters) {
    EuropeanOption::v0 = parameters[0];
    EuropeanOption::theta = parameters[1];
    EuropeanOption::rho = parameters[2];
    EuropeanOption::kappa = parameters[3];
    EuropeanOption::sigma = parameters[4];
}

unsigned Calibration::size() const {
    return SPX_options.size() + VIX_options.size();
}

std::vector<double> Calibration::SPX_Prices() const {
    std::vector<double> prices(SPX_options.size());
    for (int i = 0; i < SPX_options.size(); ++i) {
        prices[i] = SPX_options[i].Price();
    }
    return prices;
}

std::vector<double> Calibration::VIX_Prices() const {
    std::vector<double> prices(VIX_options.size());
    for (int i = 0; i < VIX_options.size(); ++i) {
        prices[i] = VIX_options[i].Price();
    }
    return prices;
}

std::vector<double> Calibration::Prices() const {
    auto SPX_prices = SPX_Prices(), VIX_prices = VIX_Prices(), prices(SPX_prices);
    prices.insert(prices.end(), VIX_prices.begin(), VIX_prices.end());
    return prices;
}



void computePrices(double *parameters, double *prices, int m, int n, void *data)
{
    // unused variables (however needed for Levmar):
    (void)m; (void)n;

    auto * calibration = static_cast<Calibration *> (data);
    Calibration::setParameters(parameters);
    for (unsigned i = 0; i < calibration->SPX_options.size(); ++i)
    {
        prices[i] = calibration->SPX_options[i].Price();
    }
    for (int i = 0; i < calibration->VIX_options.size(); ++i)
    {
        prices[i + calibration->SPX_options.size()] = calibration->VIX_options[i].Price();
    }
}

void computeGradients(double *parameters, double *gradient, int m, int n, void *data)
{
    // unused variables (however needed for Levmar):
    (void)m; (void)n;

    auto * calibration = static_cast<Calibration *> (data);
    Calibration::setParameters(parameters);

    // adding gradients for SPX options:
    for (unsigned i = 0; i < calibration->SPX_options.size(); ++i)
    {
        auto currentGradient = calibration->SPX_options[i].Jacobian();
        for (unsigned j = 0; j < EuropeanOption::nParameters; ++j)
        {
            gradient[EuropeanOption::nParameters * i + j] = currentGradient[j];
        }
    }
    // adding gradients for VIX options:
    for (unsigned i = 0; i < calibration->VIX_options.size(); ++i)
    {
        auto currentGradient = calibration->VIX_options[i].Jacobian();
        for (unsigned j = 0; j < EuropeanOption::nParameters; ++j)
        {
            gradient[ EuropeanOption::nParameters
                * (calibration->SPX_options.size() + i) + j] = currentGradient[j];
        }
    }
}

void calibrate (const Calibration &calibration, const std::string &gradientType) {

    double marketPrices[calibration.size()];
    double marketParameters[] = {EuropeanOption::v0, EuropeanOption::theta, EuropeanOption::rho,
                                 EuropeanOption::kappa, EuropeanOption::sigma};
    // Market prices:
    computePrices(marketParameters,marketPrices, (int)EuropeanOption::nParameters,
                  (int)calibration.size(), (void *) &calibration);

    // algorithm parameters
    int itmax = 300;
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0]=LM_INIT_MU;
    // stopping thresholds for
    opts[1]=1E-10;       // ||J^T e||_inf
    opts[2]=1E-10;       // ||Dp||_2
    opts[3]=1E-10;       // ||e||_2
    opts[4]= LM_DIFF_DELTA; // finite difference if used
    // >>> Enter calibrating routine >>>
    double start_s = clock();
    double p[EuropeanOption::nParameters];
    p[0] = 0.2;
    p[1] = 0.2;
    p[2] = -0.6;
    p[3] = 1.2;
    p[4] = 0.3;

    std::cout << "\r-------- -------- -------- Heston Model Calibrator -------- -------- --------"<<std::endl;
    std::cout << "Parameters:" <<  "\t" << "v0"<< "\t" << "theta" << "\t"  << "rho"
        << "\t" << "kappa" << "\t" << "sigma"<<std::endl;
    std::cout << "\rInitial point:" << "\t" << std::scientific << std::setprecision(8)
              << p[0] << "\t" << p[1] << "\t" << p[2] << "\t" <<
              p[3] << "\t" << p[4] << std::endl;

    if (gradientType == "Analytical")
        dlevmar_der(computePrices, computeGradients, p, marketPrices, (int)EuropeanOption::nParameters,
                    (int)calibration.size(), itmax, opts, info, nullptr, nullptr,
                    (void *) &calibration);

    else
        dlevmar_dif(computePrices, p, marketPrices, (int)EuropeanOption::nParameters,
                    (int)calibration.size(), itmax, opts, info, nullptr, nullptr,
                    (void *) &calibration);

    double stop_s = clock();

    // set parameters back to start:
    calibration.setParameters(marketParameters);

    // info on the calibration:
    std::cout << "Optimum found:" << std::scientific << std::setprecision(8) << "\t"<< p[0]<< "\t"
              << p[1]<< "\t"<< p[2]<< "\t"<< p[3]<< "\t"<< p[4] << std::endl;
    std::cout << "Real optimum:" << "\t" << marketParameters[0]<<"\t"<< marketParameters[1]
              << "\t"<< marketParameters[2]<< "\t"<< marketParameters[3]<< "\t"<<
              marketParameters[4] << std::endl;

    if (int(info[6]) == 6) {
        std::cout << "\r Solved: stopped by small ||e||_2 = "<< info[1] << " < " << opts[3]<< std::endl;
    } else if (int(info[6]) == 1) {
        std::cout << "\r Solved: stopped by small gradient J^T e = " << info[2] << " < " << opts[1]<< std::endl;
    } else if (int(info[6]) == 2) {
        std::cout << "\r Solved: stopped by small change Dp = " << info[3] << " < " << opts[2]<< std::endl;
    } else if (int(info[6]) == 3) {
        std::cout << "\r Unsolved: stopped by it_max " << std::endl;
    } else if (int(info[6]) == 4) {
        std::cout << "\r Unsolved: singular matrix. Restart from current p with increased mu"<< std::endl;
    } else if (int(info[6]) == 5) {
        std::cout << "\r Unsolved: no further error reduction is possible. Restart with increased mu"<< std::endl;
    } else if (int(info[6]) == 7) {
        std::cout << "\r Unsolved: stopped by invalid values, user error"<< std::endl;
    }

    std::cout << "\r-------- -------- -------- Computational cost -------- -------- --------"<<std::endl;
    std::cout << "\r          Time cost: "<< double(stop_s - start_s) /CLOCKS_PER_SEC << " seconds "<<std::endl;
    std::cout << "       # iterations: " << int(info[5]) << std::endl;
    std::cout << "# price evaluations: " << int(info[7]) << std::endl;
    std::cout << " # Jac. evaluations: "<< int(info[8]) << std::endl;
    std::cout << "# of lin sys solved: " << int(info[9])<< std::endl; //The attempts to reduce error
    std::cout << "\r-------- -------- -------- Residuals -------- -------- --------"<<std::endl;
    std::cout << "\r          ||e0||_2: " << info[0] << std::endl;
    std::cout << "          ||e*||_2: " << info[1]<<std::endl;
    std::cout << "       ||J'e||_inf: " << info[2]<<std::endl;
    std::cout << "          ||Dp||_2: " << info[3]<<std::endl;

}

bool AreSame(double a, double b)
{
    return fabs(a - b) < std::numeric_limits<double>::epsilon();
}