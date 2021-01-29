function [P_i] = FFT_Heston(x,S,K_i,r,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% European Call - Heston model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% [P_i] = FFT_Heston(K_i) - could be a vector
% INPUT: K_i = strikes - could be a vector
% OUTPUT: P_i = prices  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Npow = 15;
N = 2^(Npow); % grid point
A = 600; % upper bound
eta = A/N;
lambda = 2*pi/(N*eta); 
k = -lambda * N/2  + lambda *(0:N-1); % log-strike grid  
K = S * exp(k); % strike 
v = eta*(0:N-1);
v(1)=1e-22; %correction term: could not be equal to zero (otherwise NaN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRICING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
% Fourier transform of z_k
tr_fou = trasf_fourier(r,x,T,v);
% Trapezoidal rule
w = [0.5  ones(1,N-2)  0.5];  
h = exp(1i*(0:N-1)*pi).*tr_fou.*w*eta;
P = S * real( fft(h)/pi + max(1-exp(k-r*T),0)); % prices
%time=toc

% delete too small and too big strikes
index=find( (K>0.1*S & K<3*S) );
K=K(index); P=P(index);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PLOT
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(K,P, 'r') ;
% hold on
% axis([0  2*S  0  S]) ;
% xlabel('strike') ;
% ylabel('option price') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_i = interp1(K,P,K_i,'spline');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function fii = trasf_fourier(r,x,T,v)
fii = (exp(1i*r*v*T)).*( (characteristic_func(x,T,v-1i)-1)./(1i*v.*(1+1i*v))  );

function f = characteristic_func(x,T,u)
%--- model parameters 
epsilon=x(1); % vol-of-vol
kappa=x(2); % mean reversion speed
rho=x(3); % correlation
theta=x(4); % mean 
V0=x(5);
    
alfa = -.5*(u.*u + u*1i);
beta = kappa - rho*epsilon*u*1i;
epsilon2 = epsilon * epsilon;
gamma = .5 * epsilon2;

D = sqrt(beta .* beta - 4.0 * alfa .* gamma);

bD = beta - D;
eDt = exp(- D * T);

G = bD ./ (beta + D);
B = (bD ./ epsilon2) .* ((1.0 - eDt) ./ (1.0 - G .* eDt));
psi = (G .* eDt - 1.0) ./(G - 1.0);
A = ((kappa * theta) / (epsilon2)) * (bD * T - 2.0 * log(psi));

y = A + B*V0;

f = exp(y);

