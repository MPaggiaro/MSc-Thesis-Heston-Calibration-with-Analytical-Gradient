//
// Created by Marco Paggiaro on 19/01/2021.
//

#include "VIXOption.h"
#include "Faddeeva.hh"
#include <boost/math/quadrature/gauss.hpp>
#define xGauss boost::math::quadrature::gauss<double,64>::abscissa()
#define wGauss boost::math::quadrature::gauss<double,64>::weights()
#define I std::complex<double>(0.0,1.0)

VIXOption::VIXOption(const double K, const double T, const double r, std::string  optType,
                     std::string cfType,const double N):
        EuropeanOption(K,T,r,std::move(optType),std::move(cfType),N) { }

double VIXOption::Price() const
{
    const double a = (1 - exp(- tau_bar * kappa))/kappa,
            b = theta * (tau_bar - a),
            integralPhi = IntegralPhi();
    auto integrand = [&](std::complex<double> u)-> double
    {
        const std::complex<double> U = - u * a / tau_bar, phi = CharFunc(U);
        const std::complex<double> integrand = phi * exp(-I*u*(b + integralPhi)/tau_bar) *
                (1.0 - Faddeeva::erf(K * sqrt(- I * u)/S0)) / pow(-I * u, 3/2.0);
        return real(integrand);
    };
    double integral = 0.0;

    // Imaginary shift of u (as explained in Pacati):
    // const double zeta_T = tau_bar / ( a * pow(sigma,2)/(2*kappa) * (1 - exp(-kappa*T)));
    // with this value it explodes. Maybe we got some computation wrong.
    double zeta_T = 2;

    for (int i = 0; i < xGauss.size(); ++i)
    {
        integral += N/2 * wGauss[i] * (integrand(N/2 + N/2*xGauss[i] + zeta_T/2 * I)
                + integrand(N/2 - N/2*xGauss[i] + zeta_T/2 * I));
    }
    return S0 * exp(-r*T)/(2 * sqrt(M_PI)) * integral;
}

std::vector<double> VIXOption::Jacobian() const
{
    const double a = (1 - exp(- tau_bar * kappa))/kappa,
            b = theta * (tau_bar - a),
    // Imaginary shift of u (as explained in Pacati):
            zeta_T = 2,
            integralPhi = IntegralPhi(); // = tau_bar / ( a * pow(sigma,2)/(2*kappa));
    auto partial_integrand = [&](std::complex<double> u)-> std::complex<double>
    {
        const std::complex<double> U = - u * a / tau_bar, phi = CharFunc(U);
        const std::complex<double> integrand = S0*exp(-r*T)/(2*sqrt(M_PI)) * phi
                * exp(-I*u*(b + integralPhi) /tau_bar) *
                (1.0 - Faddeeva::erf(K * sqrt(- I * u)/S0)) / pow(-I * u, 3/2.0);
        return integrand;
    };
    std::vector<double> integral(nParameters);
    for (int i = 0; i < xGauss.size(); ++i)
    {
        auto u_up = N/2 + N/2*xGauss[i] + zeta_T/2 * I,
            u_down = N/2 - N/2*xGauss[i] + zeta_T/2 * I;
        auto phi_up = JacobianCF(u_up), phi_down = JacobianCF(u_down);
        auto part_int_up = partial_integrand(u_up), part_int_down = partial_integrand(u_down);

        for (int j = 0; j < nParameters; ++j) {
            integral[j] += N/2 * wGauss[i] * (real(phi_up[j]*part_int_up) +
                                              real(phi_down[j]*part_int_down));
        }
    }
    return integral;
}

std::complex<double> VIXOption::CharFunc(std::complex<double> u) const
{
    double sigma2 = pow(sigma,2);
    // second form:
    auto G = cosh(kappa*T/2) + (kappa - sigma2*I*u)/kappa * sinh(kappa*T/2),
            F = v0 * I * u / G * exp(-kappa*T/2);

    return pow(exp(kappa*T/2)/G, 2*kappa*theta/sigma2) * exp(F);
}

std::vector<std::complex<double>> VIXOption::JacobianCF(std::complex<double> u) const
{
    const double a = (1 - exp(- tau_bar * kappa))/kappa, sigma2 = pow(sigma,2),
            pb_pTheta = tau_bar - a,
            pa_pKappa = (tau_bar - a * (kappa * tau_bar + 1)) / kappa,
            pb_pKappa = - theta * pa_pKappa;

    auto U = [&](std::complex<double> u)-> std::complex<double>
    {
        return - u * a / tau_bar;
    };
    auto G = [&](std::complex<double> u)-> std::complex<double>
    {
        return cosh(kappa*T/2) + (kappa-sigma2*I*u)/kappa*sinh(kappa*T/2);
    };
    auto F = [&](std::complex<double> u)-> std::complex<double>
    {
        return v0 * I * u / G(u) * exp(- kappa * T / 2);
    };
    auto pG_pSigma = [&](std::complex<double> u)-> std::complex<double>
    {
        return - 2 * sigma * I * u / kappa * sinh(kappa*T/2);
    };
    auto pG_pu = - sigma2 * I / kappa * sinh(kappa*T/2);

    auto h_v0 = [&](std::complex<double> u)-> std::complex<double>
    {
        return F(u)/v0;
    };
    auto h_theta = [&](std::complex<double> u)-> std::complex<double>
    {
        return 2*kappa/sigma2 * log(exp(kappa*T/2)/G(u));
    };
    auto h_sigma = [&](std::complex<double> u)-> std::complex<double>
    {
        return - 2*theta/sigma * h_theta(u) - 2 * kappa * theta / (sigma2*G(u)) * pG_pSigma(u)
               - v0 * I * u / (pow(G(u),2)*exp(kappa*T/2)) * pG_pSigma(u);
    };
    auto h_kappa = [&](std::complex<double> u)-> std::complex<double>
    {
        return - sigma/(2*kappa) * h_sigma(u) + theta*T*I*u/(G(u)*exp(kappa*T/2))
               - v0*u*T/(2*kappa*pow(G(u),2))*(2*kappa*I + u*sigma2);
    };
    auto h_prime_kappa = [&](std::complex<double> u)-> std::complex<double>
    {
        return h_kappa(U(u)) - u / tau_bar * (F(U(u))/U(u) - 1.0/G(U(u)) * (2*kappa*theta/sigma2 + F(U(u))) * pG_pu)
                               * pa_pKappa;
    };
    std::vector<std::complex<double>> jacobian(nParameters, 0.0);

    jacobian[0] = h_v0(U(u)),
    jacobian[1] = h_theta(U(u)) - I * u / tau_bar * pb_pTheta,
    // rho derivative zero, therefore missing.
    jacobian[3] = h_prime_kappa(u) - I * u / tau_bar * pb_pKappa,
    jacobian[4] = h_sigma(U(u));

    if (nParameters > 5)
        for (int i = indexT; i < indexTPlusTau; ++i)
            jacobian[5 + i] = - I * u / tau_bar * deltaTimes[i];

    return jacobian;
}

void VIXOption::SetIndexTPlusTau()
{
    // find i such that times[i] = T + tau_bar:
    auto it = std::find_if(times.begin() + indexT, times.end(),
                           [&](double time)
                           { return fabs(T + tau_bar - time) < std::numeric_limits<double>::epsilon(); });
    indexTPlusTau = it - times.begin();
}

double VIXOption::IntegralPhi() const {
    double integralPhi = 0.0;
    for (int i = indexT; i < indexTPlusTau; ++i)
    {
        integralPhi += deltaTimes[i] * phiT[i];
    }
    return integralPhi;
}
