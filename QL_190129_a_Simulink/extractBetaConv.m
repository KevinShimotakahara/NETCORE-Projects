
%% Extract convergence of Beta
function [nSFdr1, nSFdr2, nSFdr3, nSFdr4, nSFsolar, nSFpmu,...
    dr1Data, dr2Data, dr3Data, dr4Data, solarData, pmuData] = extractBetaConv(run, beta)

global neNBs
global nD2D
global nDummy
global nUEs
global d2dScheduler

%% Convergence for one run
logfileName = strcat('./BetaResults/', num2str(neNBs), 'eNBs_', num2str(nD2D),...
    'D2Ds_', num2str(nDummy), 'Dummy_', num2str(nUEs), 'UEs_', d2dScheduler, '_Run', num2str(run),...
    '_Beta', num2str(beta), '.dat');

data = importfileLog(logfileName);

txID = data{:, 3};

% Devices' indices
drIdx = [];
for i = 1:4
   drIdx = [drIdx; find(txID == i)];
end

% Data of each device
dr1Data = table2array(data(txID == 1, 11));
dr2Data = table2array(data(txID == 2, 11));
dr3Data = table2array(data(txID == 3, 11));
dr4Data = table2array(data(txID == 4, 11));
solarData = table2array(data(txID == 7, 11));
pmuData = table2array(data(txID == 5, 11));

nSFdr1 = table2array(data(txID == 1, 9));
nSFdr2 = table2array(data(txID == 2, 9));
nSFdr3 = table2array(data(txID == 3, 9));
nSFdr4 = table2array(data(txID == 4, 9));
nSFsolar = table2array(data(txID == 7, 9));
nSFpmu = table2array(data(txID == 5, 9));

end


