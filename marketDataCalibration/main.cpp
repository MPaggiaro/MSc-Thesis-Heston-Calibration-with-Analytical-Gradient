#include <iostream>
#include <vector>
#include "EuropeanOption.h"
#include "Calibration.h"
#include <fstream>

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
                             VIX_maturities, parametersWithDisplacement, 5,
                             r, S0, false);

    double initialGuess[5 + displacement_vector.size()];
    initialGuess[0] = 0.2;
    initialGuess[1] = 0.2;
    initialGuess[2] = -0.6;
    initialGuess[3] = 1.2;
    initialGuess[4] = 0.3;
    if (!displacement_vector.empty())
        for (int i = 0; i < displacement_vector.size(); ++i)
            initialGuess[i+5] = 2e-2;

    calibrate(calibration, initialGuess);

    std::vector<std::vector<double>> convergenceParameters;
    // first parameter:
    std::vector<double> parameters (initialGuess, initialGuess + 5);
    // convergenceParameters.push_back(std::vector<double> (initialGuess, initialGuess + 5 + displacement_vector.size()));
    // std::vector<double> convergenceF;
    // save data for plots:
    for (int i = 1; i < 50; ++i) {
        calibrate(calibration, initialGuess, i);
        if (i == 1)
        {
            // add info on objective function:
            parameters.push_back(calibration.info[0]);
            convergenceParameters.push_back(parameters);
        }
        if (calibration.info[6] == 6)
            break;
        parameters = calibration.searchParameters;
        parameters.push_back(calibration.info[1]);
        convergenceParameters.push_back(parameters);
    }
    convergenceParameters.push_back(std::vector<double> (marketParameters, marketParameters+5));

    std::ofstream out("convergenceH.csv");

    for (auto& row : convergenceParameters) {
        for (auto col : row)
            out << col <<',';
        out << '\n';
    }

    Calibration calibrationDisp (SPX_strikes, SPX_maturities, VIX_strikes,
                             VIX_maturities, parametersWithDisplacement, 5 + displacement_vector.size(),
                             r, S0, true);

    calibrate(calibrationDisp, initialGuess);

    std::vector<std::vector<double>> convergenceParametersDisp;
    // first parameter:
    parameters = std::vector<double>  (initialGuess, initialGuess + 5 + displacement_vector.size());
    // convergenceParameters.push_back(std::vector<double> (initialGuess, initialGuess + 5 + displacement_vector.size()));
    // std::vector<double> convergenceF;
    // save data for plots:
    for (int i = 1; i < 50; ++i) {
        calibrate(calibrationDisp, initialGuess, i);
        if (i == 1)
        {
            // add info on objective function:
            parameters.push_back(calibrationDisp.info[0]);
            convergenceParametersDisp.push_back(parameters);
        }
        if (calibrationDisp.info[6] == 6)
            break;
        parameters = calibrationDisp.searchParameters;
        parameters.push_back(calibrationDisp.info[1]);
        convergenceParametersDisp.push_back(parameters);
    }
    convergenceParametersDisp.push_back(std::vector<double> (parametersWithDisplacement,
                                                             parametersWithDisplacement + 5 +
                                                             displacement_vector.size()));

    // add info about delta times:
    parameters = std::vector<double> (5, 0);
    parameters.insert(parameters.end(), calibrationDisp.deltaTimes.begin(),calibrationDisp.deltaTimes.end());
    convergenceParametersDisp.push_back(parameters);

    std::ofstream outHpp("convergenceH++.csv");
    for (auto& row : convergenceParametersDisp) {
        for (auto col : row)
            outHpp << col <<',';
        outHpp << '\n';
    }
    // testModel(calibration, 10);
//    testModel(noDisplacement, 10, "Numerical");

//    std::vector<double>  // CMTDates = {0,0.083333333,0.166666667,0.25,0.5,1,2},
//            // CMTRates = {0.01,	0.01, 0.02,	0.02,	0.06,	0.07,	0.15},
//            SPXStrikesMkt, SPXExpiriesMkt, SPXPricesMkt, SPXRates,
//            VIXStrikesMkt, VIXExpiriesMkt, VIXPricesMkt, VIXRates;
//
//    // get data from files:
//    std::ifstream inputSPX ("SPXData.txt"),
//        inputVIX("VIXData.txt"), inputSPXRates("ratesSPX.txt"), inputVIXRates("ratesVIX.txt");
//
//    std::string line;
//    while (std::getline(inputSPX, line))
//    {
//        std::istringstream iss(line);
//        double a, b, c;
//        if (!(iss >> a >> b >> c)) { break; }
////        std::vector<double> parts;
////        split(line, '\t', parts);
//        SPXExpiriesMkt.push_back(a);
//        SPXStrikesMkt.push_back(b);
//        SPXPricesMkt.push_back(c);
//    }
//
//    while (std::getline(inputVIX, line))
//    {
//        std::istringstream iss(line);
//        double a, b, c;
//        if (!(iss >> a >> b >> c)) { break; }
//        VIXExpiriesMkt.push_back(a);
//        VIXStrikesMkt.push_back(b);
//        VIXPricesMkt.push_back(c);
//    }
//    std::vector<std::string> SPXOptionType(SPXStrikesMkt.size(), "C");
//    for (int i = 0; i < SPXOptionType.size(); ++i)
//    {
//        if (i % 5 == 0 or i % 5 == 1)
//            SPXOptionType[i] = "P";
//    }
//
//    // These ones have to be updated "manually".
//    double S0Mkt = 3968.94;
//
//    double val;
//    while (inputSPXRates >> val)
//        SPXRates.push_back(val);
//
//    while (inputVIXRates >> val)
//        VIXRates.push_back(val);
//
//    for (int i = 0; i < SPXRates.size(); ++i) {
//        if (i % 5 == 0)
//            SPXRates[i] = SPXRates[i+2];
//        if (i % 5 == 1)
//            SPXRates[i] = SPXRates[i+1];
//    }
//
//    Calibration marketDisp(SPXStrikesMkt, SPXExpiriesMkt, SPXPricesMkt, SPXRates,
//                            VIXStrikesMkt, VIXExpiriesMkt, VIXPricesMkt, VIXRates, SPXOptionType,
//                                  S0Mkt, 100, true);
//
//    calibrateMarketData(marketDisp, initialGuess);
//    std::vector<double> bestParams = marketDisp.searchParameters;
//    double currentMin = marketDisp.info[1];
//
//    for (int i = 0; i < 10; ++i) {
//        generateGuess(initialGuess, marketDisp.nParameters);
//        calibrateMarketData(marketDisp, initialGuess);
//        if (marketDisp.info[1] < currentMin){
//            bestParams = marketDisp.searchParameters;
//            currentMin = marketDisp.info[1];
//        }
//    }
//    // save result:
//    double arrayBestParameter[bestParams.size()];
//    for (int i = 0; i < bestParams.size(); ++i) {
//        arrayBestParameter[i] = bestParams[i];
//    }
//
//    marketDisp.setParameters(arrayBestParameter);
//    auto prices = marketDisp.Prices();
//    std::ofstream outPricesHpp("heston++Prices.txt");
//    for (const auto &e : prices)
//    {
//        outPricesHpp << e << "\n";
//    }
//    std::ofstream outParamsHpp("heston++Params.txt");
//    for (const auto &e : bestParams) {
//        outParamsHpp << e << "\n";
//    }
//
//    Calibration marketNoDisp(SPXStrikesMkt, SPXExpiriesMkt, SPXPricesMkt, SPXRates,
//                           VIXStrikesMkt, VIXExpiriesMkt, VIXPricesMkt, VIXRates, SPXOptionType,
//                           S0Mkt, 100, false);
//
//    calibrateMarketData(marketNoDisp, initialGuess);
//    bestParams = marketNoDisp.searchParameters;
//    currentMin = marketNoDisp.info[1];
//
//    for (int i = 0; i < 10; ++i) {
//        generateGuess(initialGuess, marketNoDisp.nParameters);
//        calibrateMarketData(marketNoDisp, initialGuess);
//        if (marketNoDisp.info[1] < currentMin){
//            bestParams = marketNoDisp.searchParameters;
//            currentMin = marketNoDisp.info[1];
//        }
//    }
//    // save result:
//    double array2[bestParams.size()];
//    for (int i = 0; i < bestParams.size(); ++i) {
//        array2[i] = bestParams[i];
//    }
//
//    marketNoDisp.setParameters(array2);
//    prices = marketNoDisp.Prices();
//    std::ofstream outPricesH("hestonPrices.txt");
//    for (const auto &e : prices)
//    {
//        outPricesH << e << "\n";
//    }
//    std::ofstream outParamsH("hestonParams.txt");
//    for (const auto &e : bestParams) {
//        outParamsH << e << "\n";
//    }

    int a = 1;
    return 0;
}