function price = HestonPriceCarrMadan(PutCall,kappa,theta,sigma,...
    rho,v0,S0,K,T,r,q,M,N,CF,intMethod)
% Carr-Madan form
% M: upper bound; N: number of points in the grid.
dx = M/N;

if strcmp(intMethod,'TR')
    x = dx*(0:N-1); % N values.
    x(1) = 1e-22;
    % weights: trapezoidal. Maybe I can use the function trapz??
    w = ones(1,N);
    w(1) = 0.5; w(end) = w(1);
    
    int = dx*w.*CarrMadanIntegrand(x,kappa,theta,...
    sigma,rho,v0,S0,K,T,r,q,CF);
else
    % GL weights (to be fixed).
    [x,w]=lgwt(N,0,M);
    
    int = w.*CarrMadanIntegrand(x,kappa,theta,...
    sigma,rho,v0,S0,K,T,r,q,CF);
end

% Formula without alpha (as in Tankov-Cont)
price = sum(int)/pi + max(1-exp(log(K)-r*T),0);

if strcmp(PutCall,'put')
    % put price: put-call parity.
    price = price - S0*exp(-q*T) + K*exp(-r*T);
end
