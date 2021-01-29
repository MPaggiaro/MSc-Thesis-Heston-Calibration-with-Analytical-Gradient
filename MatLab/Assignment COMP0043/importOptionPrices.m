function data = importOptionPrices()
% Reads data from excel
data = readtable('COMP0043_1920_MSP_AltAsmt_Q2b.xlsx','Range','A6:K262');
data = removevars(data,{'Var6','TICKER','LastPrice','TICKER_1','LastPrice_1',...
    'STRIKE_1'});
% Changing the names of the columns (for more clarity):
data.Properties.VariableNames = {'Strike','BidCall','AskCall','BidPut','AskPut'};
end