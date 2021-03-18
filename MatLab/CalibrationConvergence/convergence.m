% Task: plot convergence results for Heston and Heston++ simulated data.

clear; close all; clc
%% 1. Get data from excel:
hData = ...
    readmatrix('convergenceH.csv');
% [daysVIX,expiryVIX, strikesVIX, bidVIX,askVIX, pricesVIX, ~,~] = ...
%     readvars('MarketData.xlsx','Sheet','Foglio3','Range','E55:L69');

hppData = readmatrix('convergenceH++.csv');

%% 3. Get data:
it = (0:size(hData,1)-2)';
itPp = (0:size(hppData,1)-3)';

v0 = abs((hData(1:end-1,1)-hData(end,1))/ hData(end,1));
v_bar = abs((hData(1:end-1,2)-hData(end,2))/ hData(end,2));
rho = abs((hData(1:end-1,3)-hData(end,3))/ hData(end,3));
kappa = abs((hData(1:end-1,4)-hData(end,4))/ hData(end,4));
sigma = abs((hData(1:end-1,5)-hData(end,5))/ hData(end,5));
r = hData(1:end-1,end-1);

v0Pp = abs((hppData(1:end-2,1)-hppData(end-1,1))/ hppData(end-1,1));
v_barPp = abs((hppData(1:end-2,2)-hppData(end-1,2))/ hppData(end-1,2));
rhoPp = abs((hppData(1:end-2,3)-hppData(end-1,3))/ hppData(end-1,3));
kappaPp = abs((hppData(1:end-2,4)-hppData(end-1,4))/ hppData(end-1,4));
sigmaPp = abs((hppData(1:end-2,5)-hppData(end-1,5))/ hppData(end-1,5));
rPp = hppData(1:end-2,end-1);

deltaPhi = zeros(size(v0Pp));
for i = 1:length(deltaPhi)
   deltaPhi(i) = sqrt(sum((hppData(end,6:end-2).*(hppData(i,6:end-2)- hppData(end-1,6:end-2))).^2)); 
end

%% 2. Plot the convergence of H:

figure('Position', [100 100 800 600]);

t = tiledlayout(1,2);
t.Padding = 'compact';
t.TileSpacing = 'compact';
title(t, 'Convergenza dell''algoritmo LM nel modello $\mathcal{H}$', 'Interpreter','latex',...
    'FontSize', 20, 'FontWeight', 'bold')

nexttile;
hold on;
xlabel ('Iterazione', 'Interpreter', 'latex')
ylabel ('Errore relativo', 'Interpreter', 'latex')
grid on
grid minor

hV0 = plot(it, v0, '-d');
hVBar = plot(it, v_bar, '-d');
hRho = plot(it, rho, '-d');
hKappa = plot(it, kappa, '-d');
hSigma = plot(it, sigma, '-d');

hLegend = legend([hV0, hVBar, hRho, hKappa, hSigma],...
    '$\left|  \frac{v_{0,i}-v_0^*}{v_0^*}  \right|$',...
    '$\left|  \frac{\overline{v}_{i}-\overline{v}^*}{\overline{v}^*}  \right|$',...
    '$\left|  \frac{\rho_i-\rho^*}{\rho^*}  \right|$',...
    '$\left|  \frac{\kappa_i - \kappa^*}{\kappa^*}  \right|$',...
    '$\left|  \frac{\sigma_i - \sigma^*}{\sigma^*}  \right|$');

set(hLegend, 'interpreter', 'latex', 'location', 'southwest', 'FontSize',15);
set([hV0, hVBar, hRho, hKappa, hSigma],'lineWidth', 1.5);

set(gca, 'YScale', 'log')
set(gca,'FontName','cmr12')

nexttile;
hR = plot(it, r, '-dk', 'lineWidth', 1.5);
hTol = yline(1e-10, '--', 'lineWidth', 1.2);
set (hTol, 'Color', [0 .5 0]);

xlabel ('Iterazione', 'Interpreter', 'latex')
ylabel ('Residuo', 'Interpreter', 'latex')
grid on
grid minor
box off
hLegend = legend([hR, hTol], '$\left\| \mathbf{r}_i \right\|$',...
    '$\varepsilon_1 = 10^{-10}$');
set(hLegend, 'interpreter', 'latex', 'location', 'southwest', 'FontSize',15);

set(gca, 'YScale', 'log')
set(gca,'FontName','cmr12')
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf,'ConvH.pdf','ContentType','vector')

%% 3: plot the convergence of H++:
figure('Position', [100 100 800 600]);

t = tiledlayout(1,2);
t.Padding = 'compact';
t.TileSpacing = 'compact';
title(t, 'Convergenza dell''algoritmo LM nel modello $\mathcal{H}++$', 'Interpreter','latex',...
    'FontSize', 20)

nexttile;
hold on;
xlabel ('Iterazione', 'Interpreter', 'latex')
ylabel ('Errore relativo', 'Interpreter', 'latex')
grid on
grid minor

hV0 = plot(it, v0, '-d');
hVBar = plot(it, v_bar, '-d');
hRho = plot(it, rho, '-d');
hKappa = plot(it, kappa, '-d');
hSigma = plot(it, sigma, '-d');

hLegend = legend([hV0, hVBar, hRho, hKappa, hSigma],...
    '$\left|  \frac{v_{0,i}-v_0^*}{v_0^*}  \right|$',...
    '$\left|  \frac{\overline{v}_{i}-\overline{v}^*}{\overline{v}^*}  \right|$',...
    '$\left|  \frac{\rho_i-\rho^*}{\rho^*}  \right|$',...
    '$\left|  \frac{\kappa_i - \kappa^*}{\kappa^*}  \right|$',...
    '$\left|  \frac{\sigma_i - \sigma^*}{\sigma^*}  \right|$');

set(hLegend, 'interpreter', 'latex', 'location', 'southwest', 'FontSize',15);
set([hV0, hVBar, hRho, hKappa, hSigma],'lineWidth', 1.5);

set(gca, 'YScale', 'log')
nexttile;
hR = plot(it, r, '-dk', 'lineWidth', 1.5);
hTol = yline(1e-10, '--', 'lineWidth', 1.2);
set (hTol, 'Color', [0 .5 0]);

xlabel ('Iterazione', 'Interpreter', 'latex')
ylabel ('Residuo', 'Interpreter', 'latex')
grid on
grid minor
hLegend = legend([hR, hTol], '$\left\| \mathbf{r}_i \right\|$',...
    '$\varepsilon_1 = 10^{-10}$');
set(hLegend, 'interpreter', 'latex', 'location', 'southwest', 'FontSize',15);

set(gca, 'YScale', 'log')
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf,'ConvH++.pdf','ContentType','vector')




