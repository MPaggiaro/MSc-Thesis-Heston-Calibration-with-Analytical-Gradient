%% Question 2B.

clear; close all; clc;

%% Import US Treasury Rates:
CMT = readtable('U.S_Treasury_yield_rates.xlsx');

% Save the settlement dates of CMT:
CMTdates = CMT.(1);
CMT = table2array(CMT(:,2:end));

CMTexpiries = [1/12 2/12 3/12 1/2 1:3 5 7 10 20 30];

% Import market data:
mktData = importOptionPrices();

settlement = datetime('06/05/2020','InputFormat','dd/MM/uuuu');
expiry = datetime('29/05/2020','InputFormat','dd/MM/uuuu');

%% Interest rate computation:
indexDate = find(CMTdates == settlement);

% 
% plot(CMTexpiries, CMT(indexDate,:));
% Spline interpolation:
T = years(expiry - settlement);
R = spline(CMTexpiries, CMT(indexDate,:),T);
tt = 0:0.1:30;
rr = spline(CMTexpiries, CMT(indexDate,:),tt);

figure
plot(CMTexpiries,CMT(indexDate,:),'o')
hold on
plot(tt,rr)
plot(T,R,'*r','LineWidth', 1.5)
hold off

%% Computation of mean price of puts and calls:
midCall = (mktData.AskCall+mktData.BidCall)/2;
midPut = (mktData.AskPut+mktData.BidPut)/2;

diffCallPut = abs(midCall-midPut);
% plot(diffCallPut)
[~,minIndex] = min(diffCallPut);
F = mktData.Strike(minIndex) + exp(R*T)*(midCall(minIndex)-midPut(minIndex));
indexK0 = find(mktData.Strike<=F,1,'last');
K0 = mktData.Strike(indexK0);

%% Selection of strikes and Q:
K = mktData.Strike;
Q = [midPut(1:indexK0-1);...
    (midPut(indexK0)+midCall(indexK0))/2;...
    midCall(indexK0+1:end)];
zeroBidPuts = find(isnan(mktData.BidPut)); 
zeroBidCalls = find(isnan(mktData.BidCall));

% In this case (manual check), we don't have strange behaviours.
% We can delete the zero bid calls and puts.
K([zeroBidPuts; zeroBidCalls]) = [];
Q([zeroBidPuts;zeroBidCalls]) = [];
    
% computation of deltaK:
deltaK = (K(3:end)-K(1:end-2))/2;
deltaK = [K(2)-K(1); deltaK; K(end)-K(end-1)];

%% Computation of VIX:
sigmaSquared = 1/T*(2*exp(R*T)*sum(deltaK.*Q./K.^2)-(F/K0-1)^2);
VIX = 100*sqrt(sigmaSquared);
% Result: 30.6493. Similar but slightly different to the one in the paper.

% Guess: what if I use the "real" F?



