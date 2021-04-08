% Task: compute values of mean data for validation:
clear; close all; clc
%% 1. Read data:
hppData = readmatrix('convergenceH++.csv');

hValData = readmatrix('modelValH.csv');
hpValData = readmatrix('modelValH++.csv');

deltaTimes = hppData(end,6:end-2);

hValData(:,end)=[];
hValData(:,12)=[];
hpValData(:,end) = [];
hpValData(:,26) = [];

%% Compute values and save them:
hInfoVal = mean(hValData);

hpInfoVal = mean(hpValData(:,1:5));
diffDisp = zeros(size(hpValData,1),1);
for i = 1:length(diffDisp)
    diffDisp(i) = sqrt(sum(hpValData(i,6:19).*deltaTimes));
end
hpInfoVal = [hpInfoVal, mean(diffDisp), mean(hpValData(:,20:28))];

dlmwrite('meanInfoValH.txt', hInfoVal);
dlmwrite('meanInfoValH++.txt', hpInfoVal);