#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <boost/math/quadrature/gauss.hpp>
using boost::math::quadrature::gauss;
#define xGauss gauss<double,64>::abscissa()
#define wGauss gauss<double,64>::weights()

bool AreSame(double a, double b);


int main() {

    double tau_bar = 0.08219178;
    std::vector<double> SPX_maturities = {
            0.08219178, 0.16438356, 0.24657534, 0.32876712, 0.4109589,  0.49315068,
            0.73972603, 0.98630137, 0.08219178, 0.16438356, 0.24657534, 0.32876712,
            0.4109589,  0.49315068, 0.73972603, 0.98630137, 0.08219178, 0.16438356,
            0.24657534, 0.32876712, 0.4109589,  0.49315068, 0.73972603, 0.98630137,
            0.08219178, 0.16438356, 0.24657534, 0.32876712, 0.4109589,  0.49315068,
            0.73972603, 0.98630137, 0.08219178, 0.16438356, 0.24657534, 0.32876712,
            0.4109589,  0.49315068, 0.73972603, 0.98630137};

    std::vector<double> VIX_maturities = {
            0.01917808, 0.01917808, 0.01917808, 0.01917808, 0.01917808, 0.03835616,
            0.03835616, 0.03835616, 0.03835616, 0.03835616, 0.05753425, 0.05753425,
            0.05753425, 0.05753425, 0.05753425, 0.08219178, 0.08219178, 0.08219178,
            0.08219178, 0.08219178, 0.16438356, 0.16438356, 0.16438356, 0.16438356,
            0.16438356, 0.24657534, 0.24657534, 0.24657534, 0.24657534, 0.24657534};

    // Erasing double elements:
    SPX_maturities.push_back(0.0);
    std::sort( SPX_maturities.begin(), SPX_maturities.end() );
    SPX_maturities.erase( std::unique(SPX_maturities.begin(), SPX_maturities.end()), SPX_maturities.end());
    std::sort( VIX_maturities.begin(), VIX_maturities.end());
    VIX_maturities.erase( std::unique(VIX_maturities.begin(), VIX_maturities.end()), VIX_maturities.end());
    std::vector<double> partialT(SPX_maturities.size() + VIX_maturities.size());
    auto it = std::set_union(SPX_maturities.begin(), SPX_maturities.end(),
                        VIX_maturities.begin(), VIX_maturities.end(), partialT.begin());
    partialT.resize(it - partialT.begin());

    auto VIX_taubar(VIX_maturities);
    for (double & i : VIX_taubar) {
        i += tau_bar;
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

    std::cout << "End of program" << std::endl;
}

bool AreSame(double a, double b)
{
    return fabs(a - b) < std::numeric_limits<double>::epsilon();
}