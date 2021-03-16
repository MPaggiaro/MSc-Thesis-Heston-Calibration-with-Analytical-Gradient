#include <iostream>
#include <vector>
#include "EuropeanOption.h"
#include "Calibration.h"

template < class T >
inline std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    for (auto ii = v.begin(); ii != v.end(); ++ii)
    {
        os << " " << *ii;
    }
    os << " ]";
    return os;
}


int main() 
{
    double r = 0.02, theta = 0.1, kappa = 3.0, sigma = 0.25, rho = -0.8, v0 = 0.08, S0 = 1.0;
    double marketParameters[] = {v0, theta, rho, kappa, sigma};
    // int nParameters = 5;
    std::vector<double> SPX_strikes = {
//            0.9371, 0.8603, 0.8112, 0.7760, 0.7470, 0.7216, 0.6699,
//            0.6137, 0.9956, 0.9868, 0.9728, 0.9588, 0.9464, 0.9358,
//            0.9175, 0.9025, 1.0427, 1.0463, 1.0499, 1.0530, 1.0562,
//            1.0593, 1.0663, 1.0766, 1.2287, 1.2399, 1.2485, 1.2659,
//            1.2646, 1.2715, 1.2859, 1.3046, 1.3939, 1.4102, 1.4291,
//            1.4456, 1.4603, 1.4736, 1.5005, 1.5328};

    0.9728, 0.9629, 0.9558, 0.9501,
0.9453, 0.9412, 0.9314, 0.9240,
0.9865, 0.9820, 0.9791, 0.9769,
0.9752, 0.9738, 0.9711, 0.9697,
1.0019, 1.0038, 1.0057, 1.0076,
1.0095, 1.0115, 1.0173, 1.0231,
1.0175, 1.0260, 1.0330, 1.0393,
1.0451, 1.0505, 1.0656, 1.0794,
1.0317, 1.0464, 1.0582, 1.0686,
1.0781, 1.0870, 1.1110, 1.1328};

    std::vector<double> SPX_maturities = {
            0.08219178, 0.16438356, 0.24657534, 0.32876712, 0.4109589,  0.49315068,
            0.73972603, 0.98630137, 0.08219178, 0.16438356, 0.24657534, 0.32876712,
            0.4109589,  0.49315068, 0.73972603, 0.98630137, 0.08219178, 0.16438356,
            0.24657534, 0.32876712, 0.4109589,  0.49315068, 0.73972603, 0.98630137,
            0.08219178, 0.16438356, 0.24657534, 0.32876712, 0.4109589,  0.49315068,
            0.73972603, 0.98630137, 0.08219178, 0.16438356, 0.24657534, 0.32876712,
            0.4109589,  0.49315068, 0.73972603, 0.98630137};

    std::vector<double> VIX_strikes = {
            0.07660333466800909, 0.07823036262505452, 0.08007866879151475,
            0.08197064388101638, 0.08371167160297435, 0.07528231789011931,
            0.07755350918528185, 0.08015741494276396, 0.08284874840358378,
            0.0853482112450969,  0.07430112345712658, 0.077055734551228,
            0.08023623852982013, 0.08354801898792094, 0.08664544590808211,
            0.0732872469748735,  0.07654634796281677, 0.08033769704753065,
            0.08431683207193616, 0.08806642128491488, 0.07084911309063555,
            0.07534553012964565, 0.08067681958626018, 0.08638533974549638,
            0.0918677586016848,  0.06910110193769552, 0.07450988995253548,
            0.08101737363349071, 0.08809320258894406, 0.0949885695945464};

    std::vector<double> VIX_maturities = {
            0.01917808, 0.01917808, 0.01917808, 0.01917808, 0.01917808, 0.03835616,
            0.03835616, 0.03835616, 0.03835616, 0.03835616, 0.05753425, 0.05753425,
            0.05753425, 0.05753425, 0.05753425, 0.08219178, 0.08219178, 0.08219178,
            0.08219178, 0.08219178, 0.16438356, 0.16438356, 0.16438356, 0.16438356,
            0.16438356, 0.24657534, 0.24657534, 0.24657534, 0.24657534, 0.24657534};

    // try new displacement initialization:
    std::vector<double> displacement_vector(14, 4e-3);

    double parametersWithDisplacement[5 + displacement_vector.size()];
    for (int i = 0; i < 5; ++i)
    {
        parametersWithDisplacement[i] = marketParameters[i];
    }
    for (int i = 0; i < displacement_vector.size(); ++i)
    {
        parametersWithDisplacement[i + 5] = displacement_vector[i];
    }

    Calibration calibration (SPX_strikes, SPX_maturities, VIX_strikes,
                             VIX_maturities, parametersWithDisplacement,
                             r, S0, true);

    auto intPhi = Calibration::IntegralDisplacement();
    auto prices = calibration.Prices();
//    auto gradients = calibration.Gradients();

    double initialGuess[EuropeanOption::nParameters];
    initialGuess[0] = 0.2;
    initialGuess[1] = 0.2;
    initialGuess[2] = -0.6;
    initialGuess[3] = 1.2;
    initialGuess[4] = 0.3;
    if (EuropeanOption::nParameters > 5)
        for (int i = 5; i < EuropeanOption::nParameters; ++i)
            initialGuess[i] = 2e-2;



//    double obFun[calibration.size()], jacFun[calibration.size() * EuropeanOption::nParameters];
//    calibration.setMarketPrices();
//    double checker[calibration.size()];
//    std::fill_n(checker, calibration.size(), -1);
//    dlevmar_chkjac(objectiveFunction, gradientObjective, initialGuess, EuropeanOption::nParameters,
//                   calibration.size(), (void *) &calibration, checker);
//    objectiveFunction(initialGuess, obFun, 1,1, (void *) &calibration);
//    gradientObjective(initialGuess, jacFun, 1, 1, (void *) &calibration);

    calibrate(calibration, initialGuess);
//    calibrate(calibration, initialGuess, true);
//    calibrate(calibration, initialGuess, false,"Numerical");

//    Calibration noDisplacement(SPX_strikes, SPX_maturities, VIX_strikes,
//                               VIX_maturities, marketParameters, r, S0);

//    calibrate(noDisplacement, initialGuess);
    testModel(calibration, 10);
//    testModel(noDisplacement, 10, "Numerical");
    return 0;
}