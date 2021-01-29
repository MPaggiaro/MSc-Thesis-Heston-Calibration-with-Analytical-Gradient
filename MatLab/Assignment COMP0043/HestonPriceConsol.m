function y = HestonPriceConsol(FlagPut,kappa,theta,sigma,...
    rho,v0,S0,K,T,r,q,lambda,Trap,Lphi,Uphi,dphi)
phi = (Lphi:dphi:Uphi);
N = length(phi);
int = zeros(N,1);
for k=1:N
int(k) = HestonProbConsol(phi(k),kappa,theta,sigma,...
    rho,v0,S0,K,T,r,q,lambda,Trap);
end
I = trapz(int)*dphi;
% The call price
HestonC = (1/2)*S0*exp(-q*T) - (1/2)*K*exp(-r*T) + I/pi;
% The put price by put-call parity
HestonP = HestonC - S0*exp(-q*T) + K*exp(-r*T);

if FlagPut == 1
    y = HestonP;
else
    y = HestonC;
end