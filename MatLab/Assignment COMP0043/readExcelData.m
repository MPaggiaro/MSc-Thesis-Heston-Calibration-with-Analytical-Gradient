function data = readExcelData()
% Reads data from excel
data = readtable('COMP0043_1920_MSP_AltAsmt_Q2a.xlsx','Range','A:H');

% Changing the names of the columns (for more clarity):
data.Properties.VariableNames = {'dateVIX','VIX','dateSPX','SPX',...
    'dateVXO','VXO','dateOEX','OEX'};
end