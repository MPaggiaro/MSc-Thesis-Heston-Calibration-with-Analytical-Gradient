function y = HestonPrice(FlagPut,kappa,theta,sigma,...
    rho,v0,S0,K,T,r,q,lambda,Trap,Lphi,Uphi,dphi)
phi = (Lphi:dphi:Uphi);
N = length(phi);
int = zeros(N,2);
for k=1:N
    int(k,:) = HestonProb(phi(k),kappa,theta,...
        sigma,rho,v0,S0,K,T,r,q,lambda,Trap);
end
I1 = trapz(int(:,1))*dphi;
I2 = trapz(int(:,2))*dphi;
P1 = 1/2 + 1/pi*I1;
P2 = 1/2 + 1/pi*I2;
HestonC = S0*exp(-q*T)*P1 - K*exp(-r*T)*P2;
HestonP = HestonC - S0*exp(-q*T) + K*exp(-r*T);

if FlagPut == 1
    y = HestonP;
else
    y = HestonC;
end