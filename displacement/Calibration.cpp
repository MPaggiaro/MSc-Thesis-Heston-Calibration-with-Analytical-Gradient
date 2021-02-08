//
// Created by Marco Paggiaro on 11/01/2021.
//

#include <iostream>
#include <iomanip>
#include <levmar/levmar.h>
#include <algorithm>
#include <iterator>
#include "Calibration.h"
#include "doublecomparison.h"
//#include <boost/math/quadrature/gauss.hpp>

Calibration::Calibration(const std::vector<double> &strikes,
                         const std::vector<double> &maturities,
                         double *parameters, double r, double S0,
                         double q, unsigned nParameters)
{
    // First, we set the parameters of the market:
    EuropeanOption::r = r;
    EuropeanOption::S0 = S0;
    EuropeanOption::q = q;
    EuropeanOption::nParameters = nParameters;
    setParameters(parameters);
    // Then, we save the information of the options:
    SPX_options.reserve(strikes.size());
    for (unsigned i = 0; i < strikes.size(); ++i)
        SPX_options.emplace_back(strikes[i], maturities[i]);
}

Calibration::Calibration(const std::vector<double> &SPX_strikes, const std::vector<double> &SPX_maturities,
                         const std::vector<double> &VIX_strikes, const std::vector<double> &VIX_maturities,
                         double *parameters, double r, double S0, bool displacementFlag, double q,
                         unsigned int nParameters) :
        Calibration(SPX_strikes, SPX_maturities, parameters, r, S0, q, nParameters)
{
    VIX_options.reserve(VIX_strikes.size());
    for (unsigned i = 0; i < VIX_strikes.size(); ++i)
        VIX_options.emplace_back(VIX_strikes[i], VIX_maturities[i]);

    if (displacementFlag)
    {
        auto SPX_t(SPX_maturities), VIX_t(VIX_maturities);
        SPX_t.push_back(0.0);
        std::sort(SPX_t.begin(), SPX_t.end());
        SPX_t.erase(std::unique(SPX_t.begin(), SPX_t.end()), SPX_t.end());
        std::sort(VIX_t.begin(), VIX_t.end());
        VIX_t.erase(std::unique(VIX_t.begin(), VIX_t.end()), VIX_t.end());
        std::vector<double> partialT(SPX_t.size() + VIX_t.size());
        auto it = std::set_union(SPX_t.begin(), SPX_t.end(),
                                 VIX_t.begin(), VIX_t.end(), partialT.begin());
        partialT.resize(it - partialT.begin());

        auto VIX_tauBar(VIX_t);
        for (double &i : VIX_tauBar)
        {
            i += EuropeanOption::tau_bar;
        }

        std::vector<double> finalT(partialT.size() + VIX_tauBar.size());
        // unite also B + tau_bar:
        auto it2 = std::set_union(partialT.begin(), partialT.end(), VIX_tauBar.begin(),
                                  VIX_tauBar.end(), finalT.begin());
        finalT.resize(it2 - finalT.begin());

        // We need to implement a function that compares two very similar doubles:
        std::vector<int> areSame;
        for (int i = 0; i < finalT.size() - 1; ++i)
        {
            if (AreSame(finalT[i], finalT[i + 1]))
                areSame.push_back(i);
        }
        for (int i : areSame)
        {
            finalT.erase(finalT.begin() + i + 1);
        }

        // compute differences:
        std::vector<double> deltaT(finalT.size() - 1);
        for (int i = 0; i < deltaT.size(); ++i)
        {
            deltaT[i] = finalT[i + 1] - finalT[i];
        }
        EuropeanOption::times = finalT;
        EuropeanOption::deltaTimes = deltaT;
        // update of the number of parameters:
        EuropeanOption::nParameters += EuropeanOption::deltaTimes.size();
        EuropeanOption::phiT.resize(EuropeanOption::deltaTimes.size());
        for (int i = 5; i < EuropeanOption::nParameters; ++i)
        {
            EuropeanOption::phiT[i - 5] = parameters[i];
            // Watch out: we haven't placed (until now) a check on the size of parameters.
            // You have to place a correct number of parameters.
        }
        // Now we are ready to set the indices for the computation of the integrals:
        for (auto & SPX_option : SPX_options)
        {
            SPX_option.SetIndexT();
        }

        for (auto & VIX_option : VIX_options)
        {
            VIX_option.SetIndexT();
            VIX_option.SetIndexTPlusTau();
        }
    }
}

void Calibration::setParameters(const double *parameters)
{
    EuropeanOption::v0 = parameters[0];
    EuropeanOption::theta = parameters[1];
    EuropeanOption::rho = parameters[2];
    EuropeanOption::kappa = parameters[3];
    EuropeanOption::sigma = parameters[4];

    if (EuropeanOption::nParameters > 5)
        // If we have to specify also the displacement, then we save it from here:
        for (int i = 5; i < EuropeanOption::nParameters; ++i)
            EuropeanOption::phiT[i - 5] = parameters[i];
}

unsigned Calibration::size() const
{
    return SPX_options.size() + VIX_options.size();
}

std::vector<double> Calibration::SPX_Prices() const
{
    std::vector<double> prices(SPX_options.size());
    for (int i = 0; i < SPX_options.size(); ++i)
    {
        prices[i] = SPX_options[i].Price();
    }
    return prices;
}

std::vector<double> Calibration::VIX_Prices() const
{
    std::vector<double> prices(VIX_options.size());
    for (int i = 0; i < VIX_options.size(); ++i)
    {
        prices[i] = VIX_options[i].Price();
    }
    return prices;
}

std::vector<double> Calibration::Prices() const
{
    auto prices = SPX_Prices(), VIX_prices = VIX_Prices();
    prices.insert(prices.end(), VIX_prices.begin(), VIX_prices.end());
    return prices;
}

std::vector<double> Calibration::SPX_Gradients() const
{
    std::vector<double> gradients;
    gradients.reserve( SPX_options.size() * EuropeanOption::nParameters);
    for (const auto & SPX_option : SPX_options)
    {
        auto currentGradient = SPX_option.Jacobian();
        gradients.insert(gradients.end(), currentGradient.begin(), currentGradient.end());
    }
    return gradients;
}

std::vector<double> Calibration::VIX_Gradients() const
{
    std::vector<double> gradients;
    gradients.reserve( VIX_options.size() * EuropeanOption::nParameters);
    for (const auto & VIX_option : VIX_options)
    {
        auto currentGradient = VIX_option.Jacobian();
        gradients.insert(gradients.end(), currentGradient.begin(), currentGradient.end());
    }
    return gradients;
}

std::vector<double> Calibration::Gradients() const
{
    auto gradients = SPX_Gradients(), VIX_gradients = VIX_Gradients();
    gradients.insert(gradients.end(), VIX_gradients.begin(), VIX_gradients.end());
    return gradients;
}

void computePrices(double *parameters, double *prices, int m, int n, void *data)
{
    // unused variables (however needed for levmar):
    (void) m;
    (void) n;

    auto *calibration = static_cast<Calibration *> (data);
    // update market parameters:
    Calibration::setParameters(parameters);
    // compute prices:
    auto pricesVector = calibration->Prices();
    // save them in the array:
    for (int i = 0; i < pricesVector.size(); ++i)
    {
        prices[i] = pricesVector[i];
    }
}

void computeGradients(double *parameters, double *gradient, int m, int n, void *data)
{
    // unused variables (however needed for Levmar):
    (void) m;
    (void) n;

    auto *calibration = static_cast<Calibration *> (data);
    // update market parameters:
    Calibration::setParameters(parameters);
    // compute gradients:
    auto gradientsVector = calibration->Gradients();
    // place them inside the array:
    for (int i = 0; i < gradientsVector.size(); ++i)
    {
        gradient[i] = gradientsVector[i];
    }
}

void calibrate(const Calibration &calibration, const double *initialGuess,
               const std::string &gradientType)
{
    // Save the "original" market parameters:
    double marketParameters[EuropeanOption::nParameters];
    marketParameters[0] = EuropeanOption::v0,
    marketParameters[1] = EuropeanOption::theta,
    marketParameters[2] = EuropeanOption::rho,
    marketParameters[3] = EuropeanOption::kappa,
    marketParameters[4] = EuropeanOption::sigma;
    if (EuropeanOption::nParameters > 5)
        for (int i = 5; i < EuropeanOption::nParameters; ++i)
            marketParameters[i] = EuropeanOption::phiT[ i - 5 ];

    // Compute market prices:
    double marketPrices[calibration.size()];
    computePrices(marketParameters, marketPrices, (int) EuropeanOption::nParameters,
                  (int) calibration.size(), (void *) &calibration);

    // algorithm parameters:
    int itMax = 300;
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0] = LM_INIT_MU;
    // stopping thresholds for
    opts[1] = 1E-10;       // ||J^T e||_inf
    opts[2] = 1E-10;       // ||Dp||_2
    opts[3] = 1E-10;       // ||e||_2
    opts[4] = LM_DIFF_DELTA; // finite difference if used
    // >>> Enter calibrating routine >>>
    double start_s = clock();

    // Search parameters:
    double p[EuropeanOption::nParameters];
//    p[0] = 0.2;
//    p[1] = 0.2;
//    p[2] = -0.6;
//    p[3] = 1.2;
//    p[4] = 0.3;
//    if (EuropeanOption::nParameters > 5)
//        for (int i = 5; i < EuropeanOption::nParameters; ++i)
//            p[i] = 2e-4;

    // Set initial guess:
    for (int i = 0; i < EuropeanOption::nParameters; ++i)
    {
        p[i] = initialGuess[i];
    }

    std::cout << "\r-------- -------- -------- Heston Model Calibrator -------- -------- --------" << std::endl;
    std::cout << "Parameters:" << "\t" << "v0" << "\t" << "theta" << "\t" << "rho"
              << "\t" << "kappa" << "\t" << "sigma" << std::endl;
    std::cout << "\rInitial point:" << "\t" << std::scientific << std::setprecision(8)
              << p[0] << "\t" << p[1] << "\t" << p[2] << "\t" <<
              p[3] << "\t" << p[4] << std::endl;

    if (gradientType == "Analytical")
        dlevmar_der(computePrices, computeGradients, p, marketPrices, (int) EuropeanOption::nParameters,
                    (int) calibration.size(), itMax, opts, info, nullptr, nullptr,
                    (void *) &calibration);

    else
        dlevmar_dif(computePrices, p, marketPrices, (int) EuropeanOption::nParameters,
                    (int) calibration.size(), itMax, opts, info, nullptr, nullptr,
                    (void *) &calibration);

    double stop_s = clock();

    // reset market parameters back to start:
    calibration.setParameters(marketParameters);

    // info on the calibration:
    std::cout << "Optimum found:" << std::scientific << std::setprecision(8) << "\t" << p[0] << "\t"
              << p[1] << "\t" << p[2] << "\t" << p[3] << "\t" << p[4] << std::endl;
    std::cout << "Real optimum:" << "\t" << marketParameters[0] << "\t" << marketParameters[1]
              << "\t" << marketParameters[2] << "\t" << marketParameters[3] << "\t" <<
              marketParameters[4] << std::endl;

    if (int(info[6]) == 6) {
        std::cout << "\r Solved: stopped by small ||e||_2 = " << info[1] << " < " << opts[3] << std::endl;
    } else if (int(info[6]) == 1) {
        std::cout << "\r Solved: stopped by small gradient J^T e = " << info[2] << " < " << opts[1] << std::endl;
    } else if (int(info[6]) == 2) {
        std::cout << "\r Solved: stopped by small change Dp = " << info[3] << " < " << opts[2] << std::endl;
    } else if (int(info[6]) == 3) {
        std::cout << "\r Unsolved: stopped by it_max " << std::endl;
    } else if (int(info[6]) == 4) {
        std::cout << "\r Unsolved: singular matrix. Restart from current p with increased mu" << std::endl;
    } else if (int(info[6]) == 5) {
        std::cout << "\r Unsolved: no further error reduction is possible. Restart with increased mu" << std::endl;
    } else if (int(info[6]) == 7) {
        std::cout << "\r Unsolved: stopped by invalid values, user error" << std::endl;
    }

    std::cout << "\r-------- -------- -------- Computational cost -------- -------- --------" << std::endl;
    std::cout << "\r          Time cost: " << double(stop_s - start_s) / CLOCKS_PER_SEC << " seconds " << std::endl;
    std::cout << "       # iterations: " << int(info[5]) << std::endl;
    std::cout << "# price evaluations: " << int(info[7]) << std::endl;
    std::cout << " # Jac. evaluations: " << int(info[8]) << std::endl;
    std::cout << "# of lin sys solved: " << int(info[9]) << std::endl; //The attempts to reduce error
    std::cout << "\r-------- -------- -------- Residuals -------- -------- --------" << std::endl;
    std::cout << "\r          ||e0||_2: " << info[0] << std::endl;
    std::cout << "          ||e*||_2: " << info[1] << std::endl;
    std::cout << "       ||J'e||_inf: " << info[2] << std::endl;
    std::cout << "          ||Dp||_2: " << info[3] << std::endl;

}

