% Exercise 1:
close all; clear;

%% a) Path simulation.
% Problem parameters:
k = 3; theta = 0.1; sigma = 0.25; rho = -0.8;
v0 = 0.08; S0 = 1; r = 0.02; q = 0;

% Discretization parameters:
T = 1; %time horizon one year.
% Number of simulated paths:
Nsim = 5;

% Number of time steps for each path:
Nsteps = 1000;

% Initial setup:
dt = T/Nsteps;
t = linspace(0,T,Nsteps+1);
S = zeros(Nsim, Nsteps+1);
v = S;
S(:,1) = S0;
v(:,1) = v0;

% Mean and covariance matrix of [Zs, Zv].
mu = [0;0];
VC = [1 rho; rho 1];


%% MC simulation of the two processes:
for i = 1:Nsteps
    % simulation of Z:
    Z = mvnrnd(mu, VC, Nsim);
    
    % Update steps:
    S(:,i+1) = S(:,i).*(1 + r*dt + sqrt(dt*v(:,i)).*Z(:,1));
    v(:,i+1) = max(v(:,i),0) + k*(theta - max(v(:,i),0))*dt + ...
        sigma*sqrt(dt*max(v(:,i),0)).*Z(:,2);
end

% Graphical representation:
figure; hold on
subplot(1,2,1); plot(t,S); title('Price');
subplot(1,2,2); plot(t,v); title('Variance');
                    
%% a1) Alternative simulation: Andersen. (useful also for Monte Carlo).


%% b) Put and call pricing with Monte Carlo:
% Strike price:
K = 1.1;

%% Simulation of the paths:
Nsim = 10e4;

S = zeros(Nsim, Nsteps+1);
v = S;
S(:,1) = S0;
v(:,1) = v0;
for i = 1:Nsteps
    % simulation of Z:
    Z = mvnrnd(mu, VC, Nsim);
    
    % Update steps:
    S(:,i+1) = S(:,i).*(1 + r*dt + sqrt(dt*v(:,i)).*Z(:,1));
    v(:,i+1) = max(v(:,i),0) + k*(theta - max(v(:,i),0))*dt + ...
        sigma*sqrt(dt*max(v(:,i),0)).*Z(:,2);
end

% Check the correct hypothesis of Q risk neutral:
[Check,~,Conf_int_Check]=normfit(S(:,end)*exp(-r*T)-S0);

% Call price:
[call_price,~,CI_call]=normfit(exp(-r*T)*...
    max(S(:,end)-K,0))

% Put price:
[put_price,~,CI_put]=normfit(exp(-r*T)*...
    max(K - S(:,end),0))


%% c) Pricing with Characteristic Function:

%% First try: maybe slow, first code in Rouah 2013.
% Problem: which value of lambda should I pass to the function????
FlagPut = 0;

% Let's try with lambda 0.
lambda = 0;
Trap = 0;
Lphi = 0.00001;
Uphi = 50;
dphi = 0.001;
tic;
y = HestonPrice(FlagPut,k,theta,sigma,...
    rho,v0,S0,K,T,r,q,lambda,Trap,Lphi,Uphi,dphi);

toc;
% It works! amazing!

% Let's try with the "consolidated technique" as in Rouah and let's see if
% it saves computational time.
tic;
y1 = HestonPriceConsol(FlagPut,k,theta,sigma,...
    rho,v0,S0,K,T,r,q,lambda,Trap,Lphi,Uphi,dphi);
toc;
% This method is faster of about 20%. 0.4 seconds compared to 0.5.
% Surely the FFT is going to be even faster.


%% Carr Madan (with and without damping factor).

format long
% Parameters:
N = 32;
M = 200;

tic
yTR = HestonPriceCarrMadan('call',k,theta,sigma,...
    rho,v0,S0,K,T,r,q,M,N,'Schoutens','TR')
toc

tic
yGL = HestonPriceCarrMadan('call',k,theta,sigma,...
    rho,v0,S0,K,T,r,q,M,N,'Schoutens','GL')
toc
% The two different forms of characteristic function give the same result.
