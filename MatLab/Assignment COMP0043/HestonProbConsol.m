function y = HestonProbConsol(phi,kappa,theta,sigma,rho,v0,S0,K,tau,r,q,lambda,Trap)
x0 = log(S0);
a = kappa*theta;
b = [kappa+lambda-rho*sigma; kappa+ lambda];
u = [.5;-.5];
d = sqrt((rho*sigma*1i*phi-b).^2 - sigma^2*(2*u*1i*phi - phi^2));
g = (b - rho*sigma*1i*phi + d) ./ (b - rho*sigma*1i*phi - d);
if Trap==1
%     c = 1/g;
%     D = (b - rho*sigma*1i*phi - d)/sigma^2*(1-exp(-d*tau)/(1-g*exp(d*tau);
%     G = (1 - c*exp(-d*tau))/(1-c);
%     C = (r-q)*i*phi*tau + a/sigma^2*((b-rho*sigma*i*phi-d)...;
elseif Trap==0
    C = (r-q)*1i*phi*tau + a/sigma^2*((b - rho*sigma*1i*phi + d)*tau ...
        - 2*log((1-g.*exp(d*tau))./(1-g)));
    D = (b - rho*sigma*1i*phi + d)/sigma^2.*(1-exp(d*tau))./(1-g.*exp(d*tau));
end
f = exp(C + D*v0 + 1i*phi*x0);
% Return the real part of the integrand
y = real(exp(-1i*phi*log(K))/1i/phi*(S0*exp(-q*tau)*f(1) - K*exp(-r*tau)*f(2)));

end