function y = CarrMadanIntegrand(u,kappa,theta,...
    sigma,rho,v0,S0,K,tau,r,q,CF)

if (strcmp(CF,'Schoutens'))
    phi = HestonCFSchoutens (u-1i,kappa,theta,...
        sigma,rho,v0,S0,tau,r,q);
else
    phi = HestonCFDelBano (u-1i,kappa,theta,...
        sigma,rho,v0,S0,tau,r,q);
end

I = exp(-1i*u*log(K)).*exp(1i*u*r*tau).*(phi - 1)./(1i*u.*(1+1i*u));
y = real(I);
end

function phi = HestonCFSchoutens (u,kappa,theta,...
        sigma,rho,v0,S0,tau,r,q)
F = S0*exp((r-q)*tau);
ksi = kappa - sigma*rho*1i*u;
d = sqrt(ksi.^2+sigma^2*(u.^2+1i*u));
g1 = (ksi+d)./(ksi-d);
g2 = 1./g1;

phi = exp(1i*u*log(F/S0) + kappa*theta/sigma^2*((ksi-d)*tau ...
    - 2*log((1-g2.*exp(-d*tau))./(1-g2))) + v0/sigma^2*(ksi-d).*(1-exp(-d*tau))...
    ./(1-g2.*exp(-d*tau)));
end

function phi = HestonCFDelBano (u,kappa,theta,...
    sigma,rho,v0,S0,tau,r,q)
F = S0*exp((r-q)*tau);
ksi = kappa - sigma*rho*1i*u;
d = sqrt(ksi.^2+sigma^2*(u.^2+1i*u));
A1 = (u.^2 + 1i*u).*sinh(d*tau/2);
A2 = d/v0.*cosh(d*tau/2) + ksi/v0.*sinh(d*tau/2);
A = A1./A2;

B = d*exp(kappa*tau/2)./(v0*A2);

phi = exp(1i*u*log(F/S0) - kappa*theta*rho*tau*1i*u/sigma - A)...
    .* B.^(2*kappa*theta/sigma^2);
    
end

