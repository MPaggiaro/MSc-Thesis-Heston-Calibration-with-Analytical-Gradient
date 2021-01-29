%% Exercise 2 VIX.
clear; close all;
%%
% Data import:
data=readExcelData;

% Computation of logarithm returns of SPX:
SPXreturns = log(data.SPX(2:end)./data.SPX(1:end-1));

%% Graphical representation: negative correlation between VIX and SPX.
figure
yyaxis left
plot(datenum(data.dateSPX(2:end)),SPXreturns)
yyaxis right
plot(datenum(data.dateVIX),data.VIX)
% Window of analysis.
xlim([datenum('01/01/08'), datenum('01-jan-2011')])
datetick('x','mmm yy','keeplimits')

% Plotting the logreturns the negative correlation is not so clear.
% Better use the actual index:
figure
yyaxis left
plot(datenum(data.dateSPX),data.SPX)
yyaxis right
plot(datenum(data.dateVIX),data.VIX)
% Window of analysis.
xlim([datenum('01/01/08'), datenum('01-jan-2011')])
datetick('x','mmm yy','keeplimits')



%% Computation of historical volatility:
% First day we can use data
firstDayIndex = 25;
lastDayIndex = find(isnat(data.dateSPX),1)-1;
dateHistVol = data.dateSPX(firstDayIndex:lastDayIndex);
HistVol = zeros(size(dateHistVol));
sampleSize = 21; %Monthly volatility.
% SPXreturns = log(data.SPX(firstDayIndex-sampleSize:lastDayIndex)./...
%     data.SPX(firstDayIndex-sampleSize-1:lastDayIndex-1));

for i = 1:length(HistVol)
    HistVol(i) = sqrt(252/sampleSize*sum(SPXreturns(firstDayIndex+i-1-sampleSize ...
        :firstDayIndex+i-1).^2));
end
% Let's make it in percentage.
HistVol = 100*HistVol;
%% Graphical representation: negative correlation between VIX and SPX.
figure
plot(datenum(dateHistVol),HistVol)
hold on
plot(datenum(data.dateVIX),data.VIX)
legend('HistVol','VIX','location','best')
% Window of analysis.
xlim([datenum('01-jan-2019'), datenum('29-apr-2020')])

datetick('x','mmm yy','keeplimits')


