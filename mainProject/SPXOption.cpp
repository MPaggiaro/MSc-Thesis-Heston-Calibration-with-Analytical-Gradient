//
// Created by Marco Paggiaro on 19/01/2021.
//

#include "SPXOption.h"
#include <boost/math/quadrature/gauss.hpp>
#define xGauss boost::math::quadrature::gauss<double,64>::abscissa()
#define wGauss boost::math::quadrature::gauss<double,64>::weights()
#define I std::complex<double>(0.0,1.0)

// Functions that calculate option price and sensitivities
double SPXOption::Price() const
{
    double integral = 0.0;
    auto integrand = [&](double u)-> double {
        std::complex<double> integrand = exp(I*u*log(S0 / K)) *
                                               CharFunc(u - 0.5 * I);
        return real(integrand) / (pow(u,2) + 0.25);
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
    auto exp_part = [&](double u)->std::complex<double> {
        // We don't have the evaluation of the characteristic function here.
        return - sqrt(S0*K)/M_PI * exp(- r * T) * exp(I * u * log(S0/K));
    };
    auto u_part = [&](double u)-> double {
        return 1.0 / (pow(u,2) + 0.25);
    };
    for (int i = 0; i < xGauss.size(); ++i)
    {
        auto u_up = N/2 + N/2*xGauss[i],
                u_down = N/2 - N/2*xGauss[i];
        // evaluation of the Jacobian vector:
        auto phi_up = JacobianCF(u_up - 0.5*I), phi_down = JacobianCF(u_down - 0.5*I);
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


