//
// Created by Marco Paggiaro on 02/01/2021.
//

#include <complex>
#include <boost/math/quadrature/gauss.hpp>
#include "EuropeanOption.h"
#define I std::complex<double>(0.0,1.0)
using boost::math::quadrature::gauss;

Market EuropeanOption::market;

double & r = EuropeanOption::market.r;
double & q = EuropeanOption::market.q;
double & kappa = EuropeanOption::market.kappa;
double & theta = EuropeanOption::market.theta;
double & sigma = EuropeanOption::market.sigma;
double & rho = EuropeanOption::market.rho;
double & v0 = EuropeanOption::market.v0;
double & S0 = EuropeanOption::market.S0;

//void initMarket(double _r, double _kappa, double _theta,
//                double _sigma, double _rho, double _v0, double S_0, double _q)
//{
//    EuropeanOption::r = _r;
//    EuropeanOption::kappa = _kappa;
//    EuropeanOption::theta = _theta;
//    EuropeanOption::sigma = _sigma;
//    EuropeanOption::rho = _rho;
//    EuropeanOption::v0 = _v0;
//    EuropeanOption::S0 = S_0;
//    EuropeanOption::q = _q;
//}

void EuropeanOption::copy(const EuropeanOption& o2)
{
    K = o2.K;
    T = o2.T;
}


EuropeanOption::~EuropeanOption() = default;

EuropeanOption& EuropeanOption::operator = (const EuropeanOption& opt2)
{ // Assignment operator (deep copy)
    if (this == &opt2) return *this;
    copy (opt2);
    return *this;
}

// Functions that calculate option price and sensitivities
double EuropeanOption::Price() const
{
    auto integrand = [&](double u)-> double {
        const std::complex<double> phi = CharFunc(u-I);
        const std::complex<double> integrand = exp(-I*u*log(K))*exp(I*u*market.r*T)
                *(phi - 1.0)/(I*u*(1.0+I*u));
        return real(integrand);
    };

    double integral = gauss<double, M>::integrate(integrand,0,N);

    double price = integral/M_PI + std::max(1-exp(log(K)-market.r*T),0.0);

    // Case put:
    if (optType == "P")
        price = price - S0*exp(-q*T) + K*exp(-r*T);

    return price;
}

std::complex<double> EuropeanOption::CharFunc(std::complex<double> u) const
{
    const double F = S0 * exp((r - q) * T);
    const std::complex<double> ksi = kappa - sigma*rho*I*u;
    const std::complex<double> d = sqrt(pow(ksi,2) + pow(sigma,2) * (pow(u,2) + I*u));

    std::complex<double> phi;

    if (cfType == "Schoutens")
    {
        const std::complex<double> g1 = (ksi + d)/(ksi - d),
                g2 = 1.0/g1;

        phi = exp(I*u*log(F/S0) + kappa*theta/pow(sigma,2)*((ksi-d)*T
                - 2.0*log((1.0-g2*exp(-d*T))/(1.0-g2))) + v0/pow(sigma,2)*(ksi-d)*
                        (1.0-exp(-d*T))/(1.0-g2*exp(-d*T)));
    }
    else if (cfType == "Cui")
    {
        const std::complex<double> A1 = (pow(u,2) + I*u)*sinh(d*(T/2)),
                A2 = d/v0*cosh(d*(T/2)) + ksi/v0*sinh(d*(T/2)),
                A = A1/A2,
                B = d*exp(kappa*T/2)/(v0*A2), D = log(B);
        phi = exp(I*u*log(F/S0) - kappa*theta*rho*T*I*u/sigma - A
                + 2*kappa*theta/pow(sigma,2)*D);
    }
    else
    {
        const std::complex<double> A1 = (pow(u,2) + I*u)*sinh(d*(T/2)),
                A2 = d/v0*cosh(d*(T/2)) + ksi/v0*sinh(d*(T/2)),
                A = A1/A2,
                B = d*exp(kappa*T/2)/(v0*A2);

        phi = exp(I*u*log(F/S0) - kappa*theta*rho*T*I*u/sigma - A)
              * pow(B,(2*kappa*theta/pow(sigma,2)));
    }
    return phi;
}

std::vector<double> EuropeanOption::ComputeJacobian() const
{
    unsigned index;
    auto integrand = [&](double v)-> double {
        const std::complex<double> integrand = 1/M_PI * exp(-I * v * log(K)) * exp(I * v * r * T) *
                JacCharFunc(v - I).at(index) / (I * v * (1.0 + I * v));
        return real(integrand);
    };
    const unsigned sizeTheta = 5;
    std::vector<double> integral;
    for (index = 0; index < sizeTheta; index ++)
        integral.push_back(gauss<double,M>::integrate(integrand,0,N));
    return integral;
}


std::vector<std::complex<double>> EuropeanOption::JacCharFunc(std::complex<double> u) const
{
    const double sigma2 = pow(sigma,2);
    const std::complex<double> ksi = kappa - sigma*rho*I*u,
        d = sqrt(pow(ksi,2) + pow(sigma,2) * (pow(u,2) + I*u)),
        A1 = (pow(u,2) + I*u)*sinh(d*(T/2)),
        A2 = d/v0*cosh(d*(T/2)) + ksi/v0*sinh(d*(T/2)),
        A = A1/A2,
        B = d*exp(kappa*T/2)/(v0*A2), D = log(B),
        // notation: dx_dy = dx/dy.
        pd_prho = - ksi * sigma * I * u / d,
        pA2_prho = - sigma * I * u * (2.0 + ksi * T) / (2.0 * d * v0) * (ksi * cosh(d * (T / 2)) + d * sinh(d * (T / 2))),
        pB_prho = exp(kappa * (T / 2)) / v0 * (1.0 / A2 * pd_prho - d / pow(A2, 2) * pA2_prho),
        pA1_prho = -I * u * (pow(u, 2) + I * u) * T * ksi * sigma / (2.0 * d) * cosh(d * (T / 2)),
        pA_prho = 1.0 / A2 * pA1_prho - A / A2 * pA2_prho,
        pB_pkappa = I / (sigma * u) * pB_prho + B * (T / 2),
        pd_psigma = (rho / sigma - 1.0 / ksi) * pd_prho + sigma * pow(u, 2) / d,
        pA1_psigma = (pow(u, 2) + I * u) * (T / 2) * pd_psigma * cosh(d * (T / 2)),
        pA2_psigma = rho / sigma * pA2_prho - (2.0 + T * ksi) / (v0 * T * ksi * I * u) * pA1_prho + sigma * T * A1 / (2 * v0),
        pA_psigma = 1.0 / A2 * pA1_psigma - A / A2 * pA2_psigma;


    const std::complex<double> h1 = - A/v0,
        h2 = 2*kappa/sigma2*D - kappa*rho*T*I*u/sigma,
        h3 = - pA_prho + 2 * kappa * theta / (sigma2 * d) * (pd_prho - d / A2 * pA2_prho) - kappa * theta * T * I * u / sigma,
        h4 = 1.0 / (sigma*I*u) * pA_prho + 2 * theta / sigma2 * D + 2 * kappa * theta / (sigma2 * B) * pB_pkappa
             - theta*rho*T*I*u/sigma,
        h5 = - pA_psigma - 4 * kappa * theta / pow(sigma, 3) * D + 2 * kappa * theta / (sigma2 * d)
                                                                   * (pd_psigma - d / A2 * pA2_psigma) + kappa * theta * rho * T * I * u / sigma2;

    const auto phi = CharFunc(u);
    std::vector<std::complex<double>> jacobian {phi*h1,phi*h2,phi*h3,phi*h4,phi*h5};
    return jacobian;
}



