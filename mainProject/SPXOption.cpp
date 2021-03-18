//
// Created by Marco Paggiaro on 19/01/2021.
//

#include "SPXOption.h"
#include <boost/math/quadrature/gauss.hpp>
#define xGauss boost::math::quadrature::gauss<double,64>::abscissa()
#define wGauss boost::math::quadrature::gauss<double,64>::weights()
#define I std::complex<double>(0.0,1.0)

SPXOption::SPXOption (const double K, const double T, const double r, std::string  optType,
                      std::string cfType,const double N):
        EuropeanOption(K,T,r,std::move(optType),std::move(cfType),N) { }

// Functions that calculate option price and sensitivities
double SPXOption::Price() const
{
    double integral = 0.0;
    auto integrand = [&](double u)-> double
    {
        std::complex<double> integrand = exp(I*u*log(S0 / K)) *
                                               CharFunc(u - 0.5 * I);
        auto pow_u = pow(u,2) + 0.25;
        return real(integrand) * exp(- pow_u * IntegralPhi() ) / pow_u;
    };

    for (int i = 0; i < xGauss.size(); ++i)
        integral += N/2*(integrand(N/2 + N/2*xGauss[i]) + integrand(N/2 - N/2*xGauss[i]))
                    *wGauss[i];

    double price = S0*exp(-q*T) - sqrt(S0*K)/M_PI * exp(-r*T)*integral;
    // Case put: put-call parity.
    if (optType == "P")
        price = price - S0*exp(-q*T) + K*exp(-r*T);
    return price;
}

std::vector<double> SPXOption::Jacobian() const
{
    std::vector<double> integral(nParameters);
    auto exp_part = [&](double u)->std::complex<double>
    {
        // We don't have the evaluation of the characteristic function here.
        return - sqrt(S0*K)/M_PI * exp(- r * T) * exp(I * u * log(S0/K));
    };
    auto u_part = [&](double u)-> double
    {
        double pow_u = pow(u,2) + 0.25;
        return exp( - pow_u * IntegralPhi() ) / pow_u;
    };
    for (int i = 0; i < xGauss.size(); ++i)
    {
        auto u_up = N/2 + N/2*xGauss[i],
                u_down = N/2 - N/2*xGauss[i];
        // evaluation of the Jacobian vector:
        auto phi_up = JacobianCF(u_up), phi_down = JacobianCF(u_down);
        auto exp_part_up = exp_part(u_up), exp_part_down = exp_part(u_down);
        auto u_part_up = u_part(u_up), u_part_down = u_part(u_down);

        for (int j = 0; j < nParameters; ++j)
            integral[j] += N/2 * wGauss[i] *
                           (real(phi_up[j] * exp_part_up) * u_part_up
                            + real(phi_down[j] * exp_part_down) * u_part_down);
    }
    return integral;
}

std::complex<double> SPXOption::CharFunc(std::complex<double> u) const
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

std::vector<std::complex<double>> SPXOption::JacobianCF(std::complex<double> u) const
{
    auto uMinusHalfI = u - 0.5 * I;
    const double sigma2 = pow(sigma,2);
    const std::complex<double> ksi = kappa - sigma * rho * I * uMinusHalfI,
            d = sqrt(pow(ksi,2) + pow(sigma,2) * (pow(uMinusHalfI, 2) + I * uMinusHalfI)),
            A1 = (pow(uMinusHalfI, 2) + I * uMinusHalfI) * sinh(d * (T / 2)),
            A2 = d/v0*cosh(d*(T/2)) + ksi/v0*sinh(d*(T/2)),
            A = A1/A2,
            B = d*exp(kappa*T/2)/(v0*A2), D = log(B),
    // notation: dx_dy = dx/dy.
    pd_pRho = - ksi * sigma * I * uMinusHalfI / d,
            pA2_pRho = - sigma * I * uMinusHalfI * (2.0 + ksi * T) / (2.0 * d * v0) * (ksi * cosh(d * (T / 2)) + d * sinh(d * (T / 2))),
            pB_pRho = exp(kappa * (T / 2)) / v0 * (1.0 / A2 * pd_pRho - d / pow(A2, 2) * pA2_pRho),
            pA1_pRho = -I * uMinusHalfI * (pow(uMinusHalfI, 2) + I * uMinusHalfI) * T * ksi * sigma / (2.0 * d) * cosh(d * (T / 2)),
            pA_pRho = 1.0 / A2 * pA1_pRho - A / A2 * pA2_pRho,
            pB_pKappa = I / (sigma * uMinusHalfI) * pB_pRho + B * (T / 2),
            pd_pSigma = (rho / sigma - 1.0 / ksi) * pd_pRho + sigma * pow(uMinusHalfI, 2) / d,
            pA1_pSigma = (pow(uMinusHalfI, 2) + I * uMinusHalfI) * (T / 2) * pd_pSigma * cosh(d * (T / 2)),
            pA2_pSigma = rho / sigma * pA2_pRho - (2.0 + T * ksi) / (v0 * T * ksi * I * uMinusHalfI) * pA1_pRho + sigma * T * A1 / (2 * v0),
            pA_pSigma = 1.0 / A2 * pA1_pSigma - A / A2 * pA2_pSigma;

    std::vector<std::complex<double>> h(5);
    h[0] = - A/v0,
    h[1] = 2*kappa/sigma2*D - kappa * rho * T * I * uMinusHalfI / sigma,
    h[2] = - pA_pRho + 2 * kappa * theta / (sigma2 * d) * (pd_pRho - d / A2 * pA2_pRho) - kappa * theta * T * I * uMinusHalfI / sigma,
    h[3] = 1.0 / (sigma * I * uMinusHalfI) * pA_pRho + 2 * theta / sigma2 * D + 2 * kappa * theta / (sigma2 * B) * pB_pKappa
           - theta * rho * T * I * uMinusHalfI / sigma,
    h[4] = - pA_pSigma - 4 * kappa * theta / pow(sigma, 3) * D + 2 * kappa * theta / (sigma2 * d)
                                                                 * (pd_pSigma - d / A2 * pA2_pSigma) + kappa * theta * rho * T * I * uMinusHalfI / sigma2;
    double pow_u = real(pow(u,2.0)) + 0.25;
    const auto phi = CharFunc(uMinusHalfI);
    std::vector<std::complex<double>> jacobian(nParameters, 0.0);
    for (int i = 0; i < 5; ++i)
    {
        jacobian[i] = h[i] * phi;
    }
    if (nParameters > 5)
        for (int i = 0; i < indexT; ++i)
            jacobian[5 + i] = - pow_u * deltaTimes[i] * phi;

    return jacobian;
}

double SPXOption::IntegralPhi() const {
    double integralPhi = 0.0;
    for (int i = 0; i < indexT; ++i)
    {
        integralPhi += deltaTimes[i] * phiT[i];
    }
    return integralPhi;
}


