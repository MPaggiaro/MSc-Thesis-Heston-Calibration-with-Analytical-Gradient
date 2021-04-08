% Task: compute and plot market implied volatilities.
% Compare them with the ones available from market.

clear; close all; clc

%% 1. Get data from excel:
[SPXType,~,daysSPX,expirySPX, strikesSPX, bidSPX,askSPX, pricesSPX, implVolSPX,~] = ...
    readvars('MarketData.xlsx','Sheet','Foglio3','Range','C8:L47');
[daysVIX,expiryVIX, strikesVIX, bidVIX,askVIX, pricesVIX, implVolVIX,~] = ...
    readvars('MarketData.xlsx','Sheet','Foglio3','Range','E55:L69');

[treasDates, treasRates] = ...
    readvars('MarketData.xlsx','Sheet','Foglio3','Range','o3:p9');

VIX_0 = 22.41;
SPX_0 = 3939.34;

SPXType = string(SPXType);
boolSPXType = SPXType == "CALL";


%% 2. plot CMTs:

% ratesSPX = interp1(treasDates,treasRates, expirySPX);
% ratesVIX = interp1(treasDates,treasRates, expiryVIX);

% It gives very wrong values. Let's try another way.

colormap jet;
cmap = colormap;
figure('Position', [100 100 500 400]);
hCMT = plot (treasDates, 100*treasRates, '-*', 'MarkerSize', 3, 'Color', cmap(1,:)); 
grid on;
grid minor;
box off;
ylabel ('Tasso d''interesse (\%)', 'Interpreter', 'latex')
xlabel ('$\tau$ (anni)', 'Interpreter', 'latex')
set(gca,'FontName','cmr12')
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf,'CMT.pdf','ContentType','vector')

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
% fileID = fopen('ratesSPX.txt','w');
% fprintf(fileID,'%f\n',ratesSPX);
% fclose(fileID);
% 
% fileID = fopen('ratesVIX.txt','w');
% fprintf(fileID,'%f\n',ratesVIX);
% fclose(fileID);


%% 3. Compute implied volatilities:
% Volatility = blsimpv(S0,Strike,Rate,Time,callPrice)
volSPX = 100*blsimpv(SPX_0, strikesSPX, ratesSPX, expirySPX, pricesSPX,'class',boolSPXType);
volVIX = 100*blsimpv(VIX_0, strikesVIX, ratesVIX, expiryVIX, pricesVIX);


%% 4. Calibrate market and import calibration results from Heston and
% Heston++ models.
fileID = fopen('Heston Without Restraints/hestonPrices.txt','r');
hPrices = fscanf(fileID, '%f');
fileID = fopen('heston++Prices.txt','r');
hppPrices = fscanf(fileID, '%f');
 

%% 5. Compute model implied volatilities.
% volSPXh = 100*blsimpv(SPX_0, strikesSPX, ratesSPX, expirySPX, hPrices(1:40),'class',boolSPXType);
% volVIXh = 100*blsimpv(VIX_0, strikesVIX, ratesVIX, expiryVIX, hPrices(41:end));
% 
% volSPXhpp = 100*blsimpv(SPX_0, strikesSPX, ratesSPX, expirySPX, hppPrices(1:40),'class',boolSPXType);
% volVIXhpp = 100*blsimpv(VIX_0, strikesVIX, ratesVIX, expiryVIX, hppPrices(41:end));


%% 6. Plot of prices:
figure('Position', [100 0 700 1600]);

SPXplot = tiledlayout(4, 2);
SPXplot.Padding = 'compact';
SPXplot.TileSpacing = 'compact';

% title(SPXplot, 'Prezzi delle opzioni S\&P500', 'Interpreter', 'latex', ...
%     'FontSize', 16, 'FontWeight', 'Bold');

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
    
    if (i == 1)
       hLegend = legend([hMktP, hHesP, hHppP, hATM],...
           'Prezzi di mercato', '$\mathcal{H}$', '$\mathcal{H}++$',...
           atmLegend, 'location', 'northeastoutside','Interpreter', 'latex');
    end
    
    set (hATM, 'LineStyle', '-.', 'Color', [0 .5 0]);
    set(hMktP, 'MarkerSize', 5);
    set([hMktP, hHesP, hHppP, hATM], 'LineWidth', 0.6);
    set(gca,'FontName','cmr12')
    
end

% Saving it (as vector image):
set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'SPXmktCal.pdf','ContentType','vector')


figure('Units', 'pixels', ...
    'Position', [100 100 700 300]);

VIXplot = tiledlayout(1,3);
VIXplot.Padding = 'compact';
VIXplot.TileSpacing = 'compact';
% title(VIXplot, 'Prezzi delle opzioni VIX', 'Interpreter', 'latex', ...
%     'FontSize', 16, 'FontWeight', 'Bold');


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
    
    if (i == 2)
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
% exportgraphics(gcf,'VIXmktCal.pdf','ContentType','vector')

%% Computation of error prices:
ehSPX = abs((hPrices(1:40) - pricesSPX)./pricesSPX);
ehppSPX = abs((hppPrices(1:40) - pricesSPX)./pricesSPX);

% ehppSPXnoOutliers = ehppSPX;
% ehppSPXnoOutliers(ehppSPXnoOutliers == max(ehppSPXnoOutliers))=[];
% ehppSPXnoOutliers(ehppSPXnoOutliers == max(ehppSPXnoOutliers))=[];
% ehppSPXnoOutliers(ehppSPXnoOutliers == max(ehppSPXnoOutliers))=[];
% ehppSPXnoOutliers(ehppSPXnoOutliers == max(ehppSPXnoOutliers))=[];

rmsreHSPX = sqrt(1/length(ehSPX)*sum(ehSPX.^2));
rmsreHppSPX = sqrt(1/length(ehppSPX)*sum(ehppSPX.^2));

ehVIX = abs((hPrices(41:end) - pricesVIX)./pricesVIX);
ehppVIX = abs((hppPrices(41:end) - pricesVIX)./pricesVIX);
rmsreHVIX = sqrt(1/length(ehVIX)*sum(ehVIX.^2));
rmsreHppVIX = sqrt(1/length(ehppVIX)*sum(ehppVIX.^2));

rmsreH = sqrt(1/length(hPrices)*(sum(ehSPX.^2)+sum(ehVIX.^2)));
rmsreHpp = sqrt(1/length(hppPrices)*(sum(ehppSPX.^2)+sum(ehppVIX.^2)));

% %% 7. Plot of errors of prices - H:
% % final error:
% figure('Position', [100 100 500 400]);
% hErrSPX = plot3 (strikesSPX, daysSPX, ehSPX, 'k.', 'MarkerSize', 10); 
% grid on;
% grid minor;
% xlabel ('$K$', 'Interpreter', 'latex')
% ylabel ('$\tau$', 'Interpreter', 'latex')
% zlabel ('Errore', 'Interpreter', 'latex')
% set (gca, 'YTick', 0:100:400);
% hold on;
% for i=1:length(strikesSPX)
%    plot3([strikesSPX(i), strikesSPX(i)], [daysSPX(i), daysSPX(i)], [0, ehSPX(i)], 'k--') 
% end
% set(gca,'FontName','cmr12')
% set(gcf, 'PaperPositionMode', 'auto');
% % exportgraphics(gcf,'finalErrorSPXHeston.pdf','ContentType','vector')
% 
% figure('Position', [100 100 500 400]);
% hErrVIX = plot3 (strikesVIX, daysVIX, ehVIX, 'k.', 'MarkerSize', 10); 
% grid on;
% grid minor;
% xlabel ('$K$', 'Interpreter', 'latex')
% ylabel ('$\tau$', 'Interpreter', 'latex')
% zlabel ('Errore', 'Interpreter', 'latex')
% set (gca, 'YTick', 0:25:100);
% hold on;
% for i=1:length(strikesVIX)
%    plot3([strikesVIX(i), strikesVIX(i)], [daysVIX(i), daysVIX(i)], [0, ehVIX(i)], 'k--') 
% end
% set(gca,'FontName','cmr12')
% set(gcf, 'PaperPositionMode', 'auto');
% % exportgraphics(gcf,'finalErrorSPXHeston.pdf','ContentType','vector')
% 
% %% 8. Plot of errors of prices - H++:
% 
% figure('Position', [100 100 500 400]);
% hppErrSPX = plot3 (strikesSPX, daysSPX, ehppSPX, 'k.', 'MarkerSize', 10); 
% grid on;
% grid minor;
% xlabel ('$K$', 'Interpreter', 'latex')
% ylabel ('$\tau$', 'Interpreter', 'latex')
% zlabel ('Errore', 'Interpreter', 'latex')
% set (gca, 'YTick', 0:100:400);
% hold on;
% for i=1:length(strikesSPX)
%    plot3([strikesSPX(i), strikesSPX(i)], [daysSPX(i), daysSPX(i)], [0, ehppSPX(i)], 'k--') 
% end
% set(gca,'FontName','cmr12')
% set(gcf, 'PaperPositionMode', 'auto');
% % exportgraphics(gcf,'finalErrorSPXHeston.pdf','ContentType','vector')
% 
% figure('Position', [100 100 500 400]);
% hppErrVIX = plot3 (strikesVIX, daysVIX, ehppVIX, 'k.', 'MarkerSize', 10); 
% grid on;
% grid minor;
% xlabel ('$K$', 'Interpreter', 'latex')
% ylabel ('$\tau$', 'Interpreter', 'latex')
% zlabel ('Errore', 'Interpreter', 'latex')
% set (gca, 'YTick', 0:25:100);
% hold on;
% for i=1:length(strikesVIX)
%    plot3([strikesVIX(i), strikesVIX(i)], [daysVIX(i), daysVIX(i)], [0, ehppVIX(i)], 'k--') 
% end
% set(gca,'FontName','cmr12')
% set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'finalErrorSPXHeston.pdf','ContentType','vector')

%% 9. Plot market implied volatilities:

figure('Position', [100 0 500 800]);

SPXplot = tiledlayout(4, 2);
SPXplot.Padding = 'compact';
SPXplot.TileSpacing = 'compact';

atmLegend = sprintf('$\\textrm{SPX}_0 = %4.2f$', SPX_0);

for i = 0:7
    t = nexttile;
    
    

    hMktP = plot(strikesSPX(5*i+1:5*(i+1)), implVolSPX(5*i+1:5*(i+1)), '-*k');
%     hBid = plot(strikesSPX(5*i+1:5*(i+1)), bidSPX(5*i+1:5*(i+1)), '^k');
%     hAsk = plot(strikesSPX(5*i+1:5*(i+1)), askSPX(5*i+1:5*(i+1)), 'vk');
    hold on
    hATM = xline(SPX_0);
        xlabel ('Strike (\$)', 'Interpreter', 'latex')
    ylabel ('Volatilit\`a (\%)', 'Interpreter', 'latex')
    grid on
    grid minor
    box off
    
    titlePlot = sprintf('$\\tau$ = %d giorni', daysSPX(5*i+1));
    title(t, titlePlot, 'Interpreter', 'latex');

%     
%     if (i == 0)
%        hLegend = legend([hMktP, hHesP, hHppP, hATM],...
%            'Prezzi di mercato', '$\mathcal{H}$', '$\mathcal{H}++$',...
%            atmLegend, 'location', 'northeast','Interpreter', 'latex');
%     end
    
    set (hATM, 'LineStyle', '-.', 'Color', [0 .5 0]);
    set(hMktP, 'MarkerSize', 5);
    set([hMktP,hATM], 'LineWidth', 0.6);
    % set(gca,'FontName','cmr12')
    
end

% Saving it (as vector image):
set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'SPXmktVol.pdf','ContentType','vector')


figure('Units', 'pixels', ...
    'Position', [100 100 700 300]);

VIXplot = tiledlayout(1,3);
VIXplot.Padding = 'compact';
VIXplot.TileSpacing = 'compact';


for i = 0:2
%     figure('Units', 'pixels', ...
%     'Position', [100 100 500 375]);
%     hold on
    
    t = nexttile;
    
     titlePlot = sprintf('$\\tau$ = %d giorni', daysSPX(5*i+1));
    title(t, titlePlot, 'Interpreter', 'latex');
    xlabel ('Strike (\$)', 'Interpreter', 'latex')
    ylabel ('Volatilit\`a (\%)', 'Interpreter', 'latex')
    grid on
    grid minor
    hold on

    % Spline for smoother model data:
%     strikes = linspace(strikesVIX(5*i+1), strikesVIX(5*(i+1)), 100);
%     sPh = makima(strikesVIX(5*i+1:5*(i+1)), hPrices(5*(i+8)+1:5*(i+9)), strikes);
%     sPhpp = makima(strikesVIX(5*i+1:5*(i+1)), hppPrices(5*(i+8)+1:5*(i+9)), strikes);
%     
    hMktP = plot(strikesVIX(5*i+1:5*(i+1)), implVolVIX(5*i+1:5*(i+1)), '-*k');
%     hBid = plot(strikesVIX(5*i+1:5*(i+1)), bidVIX(5*i+1:5*(i+1)), '^k');
%     hAsk = plot(strikesVIX(5*i+1:5*(i+1)), askVIX(5*i+1:5*(i+1)), 'vk');
%     hHesP = plot(strikes, sPh, '--r');
%     hHppP = plot(strikes, sPhpp, '-b');
    
%     if (i == 0)
%        hLegend = legend([hMktP, hHesP, hHppP],...
%            'Prezzi di mercato', '$\mathcal{H}$', '$\mathcal{H}++$',...
%            'location', 'northeast','Interpreter', 'latex');
%     end
    
    set(hMktP, 'LineWidth', 0.6);
    
    set(hMktP, 'MarkerSize', 4);
%     set(hBid, 'MarkerSize', 4);
%     set(hAsk, 'MarkerSize', 4);
    
end

set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'VIXmktVol.pdf','ContentType','vector')

%% Compute and plot displacement found:
hppData = readmatrix('timesHeston++.txt');
hppParams = readmatrix('heston++Params.txt');

deltaTimes = hppData(1,1:end-2)';
times = hppData(2,1:end-1)';
phiT = hppParams(6:end);
intPhi = [0; cumsum(phiT.*deltaTimes)];

colormap jet;
cmap = colormap;
figure('Position', [100 100 600 500]);
hIntPhi = plot (365*times, intPhi, '--.', 'MarkerSize', 10, 'Color', cmap(1,:)); 
grid on;
grid minor;
box off;
ylabel ('$I^*_\phi(0,\tau$)', 'Interpreter', 'latex')
xlabel ('$\tau$ (giorni)', 'Interpreter', 'latex')
set(gca,'FontName','cmr12')
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf,'intPhiHes++.pdf','ContentType','vector')
