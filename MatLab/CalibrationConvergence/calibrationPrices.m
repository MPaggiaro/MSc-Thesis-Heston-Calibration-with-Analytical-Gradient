% Task: plot error prices for Heston and Heston++ simulated data.

clear; close all; clc
%% 1. Get data from excel:

data = readmatrix('calibPrices.csv');
data(:,end)=[];

%% 2: Save data and compute errors:
kSPX = data(1,1:40);
kVIX = data(1,41:end);
tSPX = round(365 * data(2,1:40)); 
tVIX = round(365 * data(2,41:end));

% initial errors h model:
ieh = abs((data(3,:)-data(5,:))./ data(5,:));
ifh = abs((data(4,:)-data(5,:))./ data(5,:));

iehSPX = ieh(1:40); iehVIX = ieh(41:end);
ifhSPX = ifh(1:40); ifhVIX = ifh(41:end);

iehp = abs((data(6,:)-data(8,:))./ data(8,:));
ifhp = abs((data(7,:)-data(8,:))./ data(8,:));
iehpSPX = iehp(1:40); iehpVIX = iehp(41:end);
ifhpSPX = ifhp(1:40); ifhpVIX = ifhp(41:end);

%% 3: Plot H errors - SPX:

% initial error:
figure('Position', [100 100 500 400]);
hInitErrSPX = plot3 (kSPX, tSPX, iehSPX, 'k.', 'MarkerSize', 10); 
grid on;
grid minor;
xlabel ('$K$', 'Interpreter', 'latex')
ylabel ('$\tau$', 'Interpreter', 'latex')
zlabel ('Errore', 'Interpreter', 'latex')
set (gca, 'YTick', 0:100:400)
hold on;
for i=1:length(kSPX)
   plot3([kSPX(i), kSPX(i)], [tSPX(i), tSPX(i)], [0, iehSPX(i)], 'k--') 
end
set(gca,'FontName','cmr12')
set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'initErrorSPXHeston.pdf','ContentType','vector')

% final error:
figure('Position', [100 100 500 400]);
hFinErrSPX = plot3 (kSPX, tSPX, ifhSPX, 'k.', 'MarkerSize', 10); 
grid on;
grid minor;
xlabel ('$K$', 'Interpreter', 'latex')
ylabel ('$\tau$', 'Interpreter', 'latex')
zlabel ('Errore', 'Interpreter', 'latex')
set (gca, 'YTick', 0:100:400);
hold on;
for i=1:length(kSPX)
   plot3([kSPX(i), kSPX(i)], [tSPX(i), tSPX(i)], [0, ifhSPX(i)], 'k--') 
end
set(gca,'FontName','cmr12')
set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'finalErrorSPXHeston.pdf','ContentType','vector')

%% 4: Plot H errors - VIX:

% initial error:
figure('Position', [100 100 500 400]);
hInitErrSPX = plot3 (kVIX, tVIX, iehVIX, 'k.', 'MarkerSize', 10); 
grid on;
grid minor;
xlabel ('$K$', 'Interpreter', 'latex')
ylabel ('$\tau$', 'Interpreter', 'latex')
zlabel ('Errore', 'Interpreter', 'latex')
set (gca, 'YTick', 0:25:100)
hold on;
for i=1:length(kVIX)
   plot3([kVIX(i), kVIX(i)], [tVIX(i), tVIX(i)], [0, iehVIX(i)], 'k--') 
end
set(gca,'FontName','cmr12')
set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'initErrorVIXHeston.pdf','ContentType','vector')

% final error:
figure('Position', [100 100 500 400]);
hFinErrSPX = plot3 (kVIX, tVIX, ifhVIX, 'k.', 'MarkerSize', 10); 
grid on;
grid minor;
xlabel ('$K$', 'Interpreter', 'latex')
ylabel ('$\tau$', 'Interpreter', 'latex')
zlabel ('Errore', 'Interpreter', 'latex')
set (gca, 'YTick', 0:25:100);
hold on;
for i=1:length(kVIX)
   plot3([kVIX(i), kVIX(i)], [tVIX(i), tVIX(i)], [0, ifhVIX(i)], 'k--') 
end
set(gca,'FontName','cmr12')
set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'finalErrorVIXHeston.pdf','ContentType','vector')

%% 5: Plot H++ errors - SPX:

% initial error:
figure('Position', [100 100 500 400]);
hInitErrSPX = plot3 (kSPX, tSPX, iehpSPX, 'k.', 'MarkerSize', 10); 
grid on;
grid minor;
xlabel ('$K$', 'Interpreter', 'latex')
ylabel ('$\tau$', 'Interpreter', 'latex')
zlabel ('Errore', 'Interpreter', 'latex')
set (gca, 'YTick', 0:100:400)
hold on;
for i=1:length(kSPX)
   plot3([kSPX(i), kSPX(i)], [tSPX(i), tSPX(i)], [0, iehpSPX(i)], 'k--') 
end
set(gca,'FontName','cmr12')
set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'initErrorSPXHeston++.pdf','ContentType','vector')

% final error:
figure('Position', [100 100 500 400]);
hFinErrSPX = plot3 (kSPX, tSPX, ifhpSPX, 'k.', 'MarkerSize', 10); 
grid on;
grid minor;
xlabel ('$K$', 'Interpreter', 'latex')
ylabel ('$\tau$', 'Interpreter', 'latex')
zlabel ('Errore', 'Interpreter', 'latex')
set (gca, 'YTick', 0:100:400);
hold on;
for i=1:length(kSPX)
   plot3([kSPX(i), kSPX(i)], [tSPX(i), tSPX(i)], [0, ifhpSPX(i)], 'k--') 
end
set(gca,'FontName','cmr12')
set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'finalErrorSPXHeston++.pdf','ContentType','vector')

%% 6: Plot H++ errors - VIX:

% initial error:
figure('Position', [100 100 500 400]);
hInitErrSPX = plot3 (kVIX, tVIX, iehpVIX, 'k.', 'MarkerSize', 10); 
grid on;
grid minor;
xlabel ('$K$', 'Interpreter', 'latex')
ylabel ('$\tau$', 'Interpreter', 'latex')
zlabel ('Errore', 'Interpreter', 'latex')
set (gca, 'YTick', 0:25:100)
hold on;
for i=1:length(kVIX)
   plot3([kVIX(i), kVIX(i)], [tVIX(i), tVIX(i)], [0, iehpVIX(i)], 'k--') 
end
set(gca,'FontName','cmr12')
set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'initErrorVIXHeston++.pdf','ContentType','vector')

% final error:
figure('Position', [100 100 500 400]);
hFinErrSPX = plot3 (kVIX, tVIX, ifhpVIX, 'k.', 'MarkerSize', 10); 
grid on;
grid minor;
xlabel ('$K$', 'Interpreter', 'latex')
ylabel ('$\tau$', 'Interpreter', 'latex')
zlabel ('Errore', 'Interpreter', 'latex')
set (gca, 'YTick', 0:25:100);
hold on;
for i=1:length(kVIX)
   plot3([kVIX(i), kVIX(i)], [tVIX(i), tVIX(i)], [0, ifhpVIX(i)], 'k--') 
end
set(gca,'FontName','cmr12')
set(gcf, 'PaperPositionMode', 'auto');
% exportgraphics(gcf,'finalErrorVIXHeston++.pdf','ContentType','vector')

