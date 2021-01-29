% Program for checking the correctness of functions and gradients.

clear; close all; clc;
%%
% Characteristic function (coming from Pacati):
syms u theta kappa sigma T v0 rho;

% kappa theta sigma tau v0 
nTheta = 0.1; nKappa = 3.0; nSigma = 0.25; nRho = -0.8; nV0 = 0.08;
nT = 1;

% Original form (to be checked if it's correct):
phi_original = exp(-2*kappa*theta/sigma^2 * log(1 - 1i*u*sigma^2/(2*kappa)*...
    (1 - exp(- kappa*T))) + v0 * (1i*u*exp(- kappa*T))...
    /(1.0 - 1i*u*sigma^2/(2*kappa)*(1 - exp(- kappa*T))));

% Lin form:
G(u, theta, kappa, sigma, rho, v0, T) = cosh(kappa*T/2) + ...
    (kappa - sigma^2*1i*u)/kappa * sinh(kappa*T/2);

F(u, theta, kappa, sigma, rho, v0, T) = v0 * 1i * u / ...
    G(u, theta, kappa, sigma, rho, v0, T) * exp(-kappa*T/2);

phi_lin = (exp(kappa*T/2)/G) ^ (2*kappa*theta/sigma^2) * exp(F);

%% Graphical representation:
% figure
% subplot(2,1,1)
% % fplot([real(phi_original) real(phi_lin)], [0 200])
% fplot(real(phi_original),[0 200], 'LineWidth', 2);
% hold on;
% fplot(real(phi_lin),[0 200], '--or');
% title('Real Part')
% hold off;
% 
% subplot(2,1,2)
% fplot(imag(phi_original),[0 200], 'LineWidth', 2);
% hold on;
% fplot(imag(phi_lin),[0 200], '--or');
% title('Imaginary Part')
% hold off;
% 
% % figure 
% % fplot(real(phi_original),[0 200])
% % 
% % figure 
% % fplot(real(phi_lin),[0 200])
% 
% % exp(-2*kappa*theta/sigma2 * log(1.0 - I*u*sigma2/(2*kappa)*(1 - ekt))
% % //            );
% 
% % double(real(phi_original(1)))
% % double(real(phi_lin(1)))

%% Check the derivatives:
%% v0: OK

Dphi_dv0 = diff(phi_lin, v0)

myderv0 = phi_lin * F / v0

%% theta:

Dphi_dtheta = diff(phi_lin, theta)
h_theta = 2 * kappa/sigma^2 * log(exp(kappa*T/2)/G);
mydertheta = phi_lin * h_theta
% OK

% Graphical representation:
fplot(real(Dphi_dtheta(u,nTheta, nKappa, nSigma, nRho, nV0, nT)),[0 200], 'LineWidth', 2)
hold on
fplot(real(mydertheta(u,nTheta, nKappa, nSigma, nRho, nV0, nT)),[0 200], '--or')
%% sigma:
Dphi_dsigma = diff(phi_lin, sigma)

dG_dsigma = - 2 * sigma * 1i * u / kappa * sinh(kappa*T/2);
hsigma = - 2*theta/sigma *  h_theta - ...
    2*kappa*theta/(sigma^2*G)*dG_dsigma  - v0*1i*u/(G^2*exp(kappa*T/2))*dG_dsigma;
mydersigma = phi_lin * hsigma
% Graphical representation:
close all
fplot(real(Dphi_dsigma(u,nTheta, nKappa, nSigma, nRho, nV0, nT)),[0 200],'LineWidth', 2)
hold on
fplot(real(mydersigma(u,nTheta, nKappa, nSigma, nRho, nV0, nT)),[0 200], '--or')

% OK.

%% kappa:
Dphi_dkappa = diff(phi_lin, kappa)

hkappa = - sigma/(2*kappa);

% Problem solved!
