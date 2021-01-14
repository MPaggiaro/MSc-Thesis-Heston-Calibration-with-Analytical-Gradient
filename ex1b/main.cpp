#include <iostream>
#include <ctime>
#include <vector>
#include "EuropeanOption.h"
#include "Calibration.h"
#include <levmar/levmar.h>
#include <iomanip>

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

void computePrices(double *p, double *x, int m, int n, void *data)
{
    auto * calibration = static_cast<Calibration *> (data);
    Calibration::setParameters(p);
    for (unsigned i = 0; i < calibration->options.size(); ++i)
    {
        x[i] = calibration->options[i].Price();
    }
}

void computeGradients(double *p, double *jac, int m, int n, void *data)
{
    auto * calibration = static_cast<Calibration *> (data);
    Calibration::setParameters(p);
    for (unsigned i = 0; i < n; ++i)
    {
        auto currentGradient = calibration->options[i].Jacobian();
        for (unsigned j = 0; j < m; ++j)
        {
            jac[m * i + j] = currentGradient[j];
        }
    }
}


int main() 
{
    double r = 0.02, theta = 0.1, kappa = 3.0, sigma = 0.25, rho = -0.8, v0 = 0.08, S0 = 1.0;
    double marketParameters[] = {v0, theta, rho, kappa, sigma};
    int nParameters = 5;
    std::vector<double> strikes = {
            0.9371, 0.8603, 0.8112, 0.7760, 0.7470, 0.7216, 0.6699, 0.6137,
            0.9956, 0.9868, 0.9728, 0.9588, 0.9464, 0.9358, 0.9175, 0.9025,
            1.0427, 1.0463, 1.0499, 1.0530, 1.0562, 1.0593, 1.0663, 1.0766,
            1.2287, 1.2399, 1.2485, 1.2659, 1.2646, 1.2715, 1.2859, 1.3046,
            1.3939, 1.4102, 1.4291, 1.4456, 1.4603, 1.4736, 1.5005, 1.5328};

    // currently maturities and strikes have the same size.
    std::vector<double> maturities = {
            0.119047619047619,  0.238095238095238,	0.357142857142857, 0.476190476190476,	0.595238095238095, 0.714285714285714, 1.07142857142857, 1.42857142857143,
            0.119047619047619,  0.238095238095238,  0.357142857142857, 0.476190476190476, 0.595238095238095, 0.714285714285714	,1.07142857142857, 1.42857142857143	,
            0.119047619047619, 	0.238095238095238,	0.357142857142857,	0.476190476190476,	0.595238095238095,	0.714285714285714,	1.07142857142857,	1.42857142857143,
            0.119047619047619,	0.238095238095238,	0.357142857142857,	0.476190476190476	,0.595238095238095,	0.714285714285714,	1.07142857142857,	1.42857142857143,
            0.119047619047619,	0.238095238095238,	0.357142857142857,	0.476190476190476,	0.595238095238095,	0.714285714285714,	1.07142857142857,	1.42857142857143};

    // Calibration calibration(strikes, maturities, 0.02,0,3.0,0.1,0.25,-0.8,0.08,1);
    Calibration calibration(strikes, maturities, marketParameters, r, S0);
//    double start = std::clock();
//    auto prices = calibration.GetPrices();
//    double stop = std::clock();
//    std::cout << "Elapsed time: " << double(stop - start)/CLOCKS_PER_SEC << " seconds" << std::endl;
//    start = std::clock();
//    auto gradients = calibration.GetGradients();
//    stop = std::clock();

    double marketPrices[strikes.size()];
    computePrices(marketParameters,marketPrices,nParameters,strikes.size(),
                  (void *) &calibration);

    // algorithm parameters
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0]=LM_INIT_MU;
    // stopping thresholds for
    opts[1]=1E-10;       // ||J^T e||_inf
    opts[2]=1E-10;       // ||Dp||_2
    opts[3]=1E-10;       // ||e||_2
    opts[4]= LM_DIFF_DELTA; // finite difference if used

    // >>> Enter calibrating routine >>>
    double start_s = clock();
    double p[nParameters];
    p[0] = 0.2;
    p[1] = 0.2;
    p[2] = -0.6;
    p[3] = 1.2;
    p[4] = 0.3;

    std::cout << "\r-------- -------- -------- Heston Model Calibrator -------- -------- --------"<<std::endl;
    std::cout << "Parameters:" <<  "\t     v0"<<"\t     theta"<<  "\t   <<   rho"  "\t         kappa"<< "\t       sigma"<<std::endl;
    std::cout << "\r Initial point:" << "\t" << std::scientific << std::setprecision(8)
              << p[0] << "\t" << p[1] << "\t" << p[2] << "\t" <<
              p[3] << "\t" << p[4] << std::endl;

    double jac[nParameters*strikes.size()];
    dlevmar_der(computePrices, computeGradients, p, marketPrices, nParameters,
                strikes.size(), 100, opts, info, nullptr, nullptr, (void *) &calibration);

    double stop_s = clock();
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
        std::cout << "\r Unsolved: stopped by itmax " << std::endl;
    } else if (int(info[6]) == 4) {
        std::cout << "\r Unsolved: singular matrix. Restart from current p with increased mu"<< std::endl;
    } else if (int(info[6]) == 5) {
        std::cout << "\r Unsolved: no further error reduction is possible. Restart with increased mu"<< std::endl;
    } else if (int(info[6]) == 7) {
        std::cout << "\r Unsolved: stopped by invalid values, user error"<< std::endl;
    }

    std::cout << "\r-------- -------- -------- Computational cost -------- -------- --------"<<std::endl;
    std::cout << "\r          Time cost: "<< double(stop_s - start_s) /CLOCKS_PER_SEC << " seconds "<<std::endl;
    std::cout << "         Iterations: " << int(info[5]) << std::endl;
    std::cout << "         pv  Evalue: " << int(info[7]) << std::endl;
    std::cout << "         Jac Evalue: "<< int(info[8]) << std::endl;
    std::cout << "# of lin sys solved: " << int(info[9])<< std::endl; //The attempts to reduce error
    std::cout << "\r-------- -------- -------- Residuals -------- -------- --------"<<std::endl;
    std::cout << "\r          ||e0||_2: " << info[0] << std::endl;
    std::cout << "          ||e*||_2: " << info[1]<<std::endl;
    std::cout << "       ||J'e||_inf: " << info[2]<<std::endl;
    std::cout << "          ||Dp||_2: " << info[3]<<std::endl;
    return 0;
}