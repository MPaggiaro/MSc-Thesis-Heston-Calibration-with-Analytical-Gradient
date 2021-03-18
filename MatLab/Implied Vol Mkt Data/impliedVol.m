% Task: compute and plot market implied volatilities.
% Compare them with the ones available from market.

clear; close all; clc

%% 1. Get data from excel:
[SPXType,~,daysSPX,expirySPX, strikesSPX, bidSPX,askSPX, pricesSPX, implVolSPX,~] = ...
    readvars('MarketData.xlsx','Sheet','Foglio3','Range','C8:L47');
[daysVIX,expiryVIX, strikesVIX, bidVIX,askVIX, pricesVIX, ~,~] = ...
    readvars('MarketData.xlsx','Sheet','Foglio3','Range','E55:L69');

% [treasDates, treasRates] = ...
%     readvars('MarketData.xlsx','Sheet','Foglio3','Range','o3:p9');

VIX_0 = 22.41;
SPX_0 = 3939.34;

SPXType = string(SPXType);
boolSPXType = SPXType == "CALL";


%% 2. Interpolate interest rates from CMTs:

% ratesSPX = interp1(treasDates,treasRates, expirySPX);
% ratesVIX = interp1(treasDates,treasRates, expiryVIX);

% It gives very wrong values. Let's try another way.

%% Compute market rates from implied volatilities:

% first, compute SPX rates from implied volatilities:
ratesSPX = zeros(size(implVolSPX));
for i = 1:length(ratesSPX)
    f = @(r) implVolSPX(i) - 100*blsimpv(SPX_0, strikesSPX(i), r, expirySPX(i), pricesSPX(i),'class',boolSPXType(i));
    ratesSPX(i) = fzero(f,0.01);
end

% Then, select the rates for VIX rates interpolation.
referenceRates = zeros(length(pricesSPX)/5, 1);
for i = 1: length(referenceRates)
   referenceRates(i) = mean(ratesSPX(5*(i-1) + 3: 5*(i-1)+5)); 
end

referenceDates = expirySPX(1:5:36);

ratesVIX = interp1(referenceDates, referenceRates, expiryVIX);
% Save rates (in order to be used in C++ calibration):
fileID = fopen('ratesSPX.txt','w');
fprintf(fileID,'%f\n',ratesSPX);
fclose(fileID);

fileID = fopen('ratesVIX.txt','w');
fprintf(fileID,'%f\n',ratesVIX);
fclose(fileID);


%% 3. Compute implied volatilities:
% Volatility = blsimpv(S0,Strike,Rate,Time,callPrice)
volSPX = 100*blsimpv(SPX_0, strikesSPX, ratesSPX, expirySPX, pricesSPX,'class',boolSPXType);
volVIX = 100*blsimpv(VIX_0, strikesVIX, ratesVIX, expiryVIX, pricesVIX);


%% 4. Calibrate market and import calibration results from Heston and
% Heston++ models.
fileID = fopen('hestonPrices.txt','r');
hPrices = fscanf(fileID, '%f');
fileID = fopen('heston++Prices.txt','r');
hppPrices = fscanf(fileID, '%f');
 

%% 5. Compute model implied volatilities.
volSPXh = 100*blsimpv(SPX_0, strikesSPX, ratesSPX, expirySPX, hPrices(1:40),'class',boolSPXType);
volVIXh = 100*blsimpv(VIX_0, strikesVIX, ratesVIX, expiryVIX, hPrices(41:end));

volSPXhpp = 100*blsimpv(SPX_0, strikesSPX, ratesSPX, expirySPX, hppPrices(1:40),'class',boolSPXType);
volVIXhpp = 100*blsimpv(VIX_0, strikesVIX, ratesVIX, expiryVIX, hppPrices(41:end));

% %% 5. Plot implied volatilities (as in Pacati et al.).
% 
% % We want smoother data for splines:
% for i = 0:7
%     figure('Units', 'pixels', ...
%     'Position', [100 100 500 375]);
%     hold on
%     
%     % Spline for smoother model data:
%     strikes = linspace(strikesSPX(5*i+1), strikesSPX(5*(i+1)), 100);
%     sVh = pchip(strikesSPX(5*i+1:5*(i+1)), volSPXh(5*i+1:5*(i+1)), strikes);
%     sVhpp = pchip(strikesSPX(5*i+1:5*(i+1)), volSPXhpp(5*i+1:5*(i+1)), strikes);
%     
%     hMktIV = plot(strikesSPX(5*i+1:5*(i+1)), volSPX(5*i+1:5*(i+1)), '*k');
%     hHesIV = plot(strikes, sVh, '--r');
%     hHppIV = plot(strikes, sVhpp, '-b');
%     
% end
% 
% for i = 0:2
%     figure('Units', 'pixels', ...
%     'Position', [100 100 500 375]);
%     hold on
%     
%     % Spline for smoother model data:
%     strikes = linspace(strikesVIX(5*i+1), strikesVIX(5*(i+1)), 100);
%     sVh = pchip(strikesVIX(5*i+1:5*(i+1)), volVIXh(5*i+1:5*(i+1)), strikes);
%     sVhpp = pchip(strikesVIX(5*i+1:5*(i+1)), volVIXhpp(5*i+1:5*(i+1)), strikes);
%     
%     hMktIV = plot(strikesVIX(5*i+1:5*(i+1)), volVIX(5*i+1:5*(i+1)), '*k');
%     hHesIV = plot(strikes, sVh, '--r');
%     hHppIV = plot(strikes, sVhpp, '-b');
% end

%% 6. Plot of prices:



figure('Position', [100 0 700 1600]);


SPXplot = tiledlayout(4, 2);
SPXplot.Padding = 'compact';
SPXplot.TileSpacing = 'compact';

title(SPXplot, 'Prezzi delle opzioni S\&P500', 'Interpreter', 'latex', ...
    'FontSize', 16, 'FontWeight', 'Bold');

atmLegend = sprintf('$\\textrm{SPX}_0 = %4.2f$', SPX_0);

for i = 0:7
    t = nexttile;
    
    titlePlot = sprintf('$\\tau$ = %d giorni', daysSPX(5*i+1));
    title(t, titlePlot, 'Interpreter', 'latex');
    xlabel ('Strike (\$)', 'Interpreter', 'latex')
    ylabel ('Prezzo (\$)', 'Interpreter', 'latex')
    grid on
    grid minor

    hold on
    % Spline for smoother model data:
    strikes = linspace(strikesSPX(5*i+1), strikesSPX(5*(i+1)), 100);
    sPh = pchip(strikesSPX(5*i+1:5*(i+1)), hPrices(5*i+1:5*(i+1)), strikes);
    sPhpp = pchip(strikesSPX(5*i+1:5*(i+1)), hppPrices(5*i+1:5*(i+1)), strikes);
    
    hMktP = plot(strikesSPX(5*i+1:5*(i+1)), pricesSPX(5*i+1:5*(i+1)), '*k');
%     hBid = plot(strikesSPX(5*i+1:5*(i+1)), bidSPX(5*i+1:5*(i+1)), '^k');
%     hAsk = plot(strikesSPX(5*i+1:5*(i+1)), askSPX(5*i+1:5*(i+1)), 'vk');
    hHesP = plot(strikes, sPh, '--r');
    hHppP = plot(strikes, sPhpp, '-b');
    hATM = xline(SPX_0);
    
    if (i == 0)
       hLegend = legend([hMktP, hHesP, hHppP, hATM],...
           'Prezzi di mercato', '$\mathcal{H}$', '$\mathcal{H}++$',...
           atmLegend, 'location', 'northeast','Interpreter', 'latex');
    end
    
    set (hATM, 'LineStyle', '-.', 'Color', [0 .5 0]);
    set(hMktP, 'MarkerSize', 5);
    set([hMktP, hHesP, hHppP, hATM], 'LineWidth', 0.6);
    set(gca,'FontName','cmr12')
    
end

% Saving it (as vector image):
% set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'SPXmktCal.pdf','ContentType','vector')


figure('Units', 'pixels', ...
    'Position', [100 100 700 300]);

VIXplot = tiledlayout(1,3);
VIXplot.Padding = 'compact';
VIXplot.TileSpacing = 'compact';
title(VIXplot, 'Prezzi delle opzioni VIX', 'Interpreter', 'latex', ...
    'FontSize', 16, 'FontWeight', 'Bold');


for i = 0:2
%     figure('Units', 'pixels', ...
%     'Position', [100 100 500 375]);
%     hold on
    
    t = nexttile;
    
    titlePlot = sprintf('$\\tau$ = %d giorni', daysVIX(5*i+1));
    title(t, titlePlot, 'Interpreter', 'latex');
    xlabel ('Strike (\$)', 'Interpreter', 'latex')
    ylabel ('Prezzo (\$)', 'Interpreter', 'latex')
    grid on
    grid minor
    
    hold on
    % Spline for smoother model data:
    strikes = linspace(strikesVIX(5*i+1), strikesVIX(5*(i+1)), 100);
    sPh = makima(strikesVIX(5*i+1:5*(i+1)), hPrices(5*(i+8)+1:5*(i+9)), strikes);
    sPhpp = makima(strikesVIX(5*i+1:5*(i+1)), hppPrices(5*(i+8)+1:5*(i+9)), strikes);
    
    hMktP = plot(strikesVIX(5*i+1:5*(i+1)), pricesVIX(5*i+1:5*(i+1)), '*k');
    hBid = plot(strikesVIX(5*i+1:5*(i+1)), bidVIX(5*i+1:5*(i+1)), '^k');
    hAsk = plot(strikesVIX(5*i+1:5*(i+1)), askVIX(5*i+1:5*(i+1)), 'vk');
    hHesP = plot(strikes, sPh, '--r');
    hHppP = plot(strikes, sPhpp, '-b');
    
    if (i == 0)
       hLegend = legend([hMktP, hHesP, hHppP],...
           'Prezzi di mercato', '$\mathcal{H}$', '$\mathcal{H}++$',...
           'location', 'northeast','Interpreter', 'latex');
    end
    
    set([hMktP, hHesP, hHppP], 'LineWidth', 0.6);
    
    set(hMktP, 'MarkerSize', 4);
    set(hBid, 'MarkerSize', 4);
    set(hAsk, 'MarkerSize', 4);
    
end

set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf,'VIXmktCal.pdf','ContentType','vector')

