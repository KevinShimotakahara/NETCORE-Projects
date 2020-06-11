%% Function to extract the results
function [drDelay, solarDelay, pmuDelay,...
    drTpt, solarTpt, pmuTpt,...
    drPDR, solarPDR, pmuPDR,...   
    drReward, solarReward, pmuReward,...
    drCQI, solarCQI, pmuCQI,...
    drDelayCI, solarDelayCI, pmuDelayCI,...
    drPDRCI, solarPDRCI, pmuPDRCI,...
    drCQICI, solarCQICI, pmuCQICI] =...
    extractBetaResults(nRuns, beta)



%% Global variables
global neNBs
global nD2D
global nDummy
global nUEs
global d2dScheduler

explTime = 7000;

dr1AvgDelay = zeros(nRuns, 1);
dr2AvgDelay = zeros(nRuns, 1);
dr3AvgDelay = zeros(nRuns, 1);
dr4AvgDelay = zeros(nRuns, 1);
solarAvgDelay = zeros(nRuns, 1);
pmuAvgDelay = zeros(nRuns, 1);

dr1AvgTpt = zeros(nRuns, 1);
dr2AvgTpt = zeros(nRuns, 1);
dr3AvgTpt = zeros(nRuns, 1);
dr4AvgTpt = zeros(nRuns, 1);
solarAvgTpt = zeros(nRuns, 1);
pmuAvgTpt = zeros(nRuns, 1);

dr1AvgPDR = zeros(nRuns, 1);
dr2AvgPDR = zeros(nRuns, 1);
dr3AvgPDR = zeros(nRuns, 1);
dr4AvgPDR = zeros(nRuns, 1);
solarAvgPDR = zeros(nRuns, 1);
pmuAvgPDR = zeros(nRuns, 1);

% dr2DelayCI = [];
% dr3DelayCI = [];
% dr4DelayCI = [];
% solarDelayCI = [];
% pmuDelayCI = [];
% 
% dr2TptCI = [];
% dr3TptCI = [];
% dr4TptCI = [];
% solarTptCI = [];
% pmuTptCI = [];
% *********************************************************************** %

%% Delay/Throughput/PDR computations
for r = 1:nRuns
    logfileName = strcat('./BetaResults/', num2str(neNBs), 'eNBs_', num2str(nD2D),...
        'D2Ds_', num2str(nDummy), 'Dummy_', num2str(nUEs), 'UEs_', d2dScheduler, '_Run', num2str(r),...
        '_Beta', num2str(beta), '.dat');

    data = importfileLog(logfileName);
    
    txID = data{:, 3};
    
    % Devices' indices
    drIdx = [];
    for i = 1:4
       drIdx = [drIdx; find(txID == i)];
    end
    
    % Data of each device
    dr1Data = data(txID == 1, :);
    dr2Data = data(txID == 2, :);
    dr3Data = data(txID == 3, :);
    dr4Data = data(txID == 4, :);
    solarData = data(txID == 7, :);
    pmuData = data(txID == 5, :);

    dr1Data = dr1Data((table2array(dr1Data(:, 8)) > explTime), :); 
    dr2Data = dr2Data((table2array(dr2Data(:, 8)) > explTime), :); 
    dr3Data = dr3Data((table2array(dr3Data(:, 8)) > explTime), :); 
    dr4Data = dr4Data((table2array(dr4Data(:, 8)) > explTime), :); 
    solarData = solarData((table2array(solarData(:, 8)) > explTime), :); 
    pmuData = pmuData((table2array(pmuData(:, 8)) > explTime), :); 
    
    % Delay results
%     dr1AvgDelay(r) = mean(table2array(dr1Data(:, 10)) - table2array(dr1Data(:, 8)));
    dr2AvgDelay(r) = mean(table2array(dr2Data(:, 10)) - table2array(dr2Data(:, 8)));
    dr3AvgDelay(r) = mean(table2array(dr3Data(:, 10)) - table2array(dr3Data(:, 8)));
    dr4AvgDelay(r) = mean(table2array(dr4Data(:, 10)) - table2array(dr4Data(:, 8)));
    solarAvgDelay(r) = mean(table2array(solarData(:, 10)) - table2array(solarData(:, 8)));
    pmuAvgDelay(r) = mean(table2array(pmuData(:, 10)) - table2array(pmuData(:, 8)));
    
    % Instantaneous throughput results
%     dr1AvgTpt(r) = mean(table2array(dr1Data(:, 7)) ./ (dr1AvgDelay(r) * 1e-3));
    dr2AvgTpt(r) = mean(table2array(dr2Data(:, 7)) ./ (dr2AvgDelay(r) * 1e-3));
    dr3AvgTpt(r) = mean(table2array(dr3Data(:, 7)) ./ (dr3AvgDelay(r) * 1e-3));
    dr4AvgTpt(r) = mean(table2array(dr4Data(:, 7)) ./ (dr4AvgDelay(r) * 1e-3));
    solarAvgTpt(r) = mean(table2array(solarData(:, 7)) ./ (solarAvgDelay(r) * 1e-3));
    pmuAvgTpt(r) = mean(table2array(pmuData(:, 7)) ./ (pmuAvgDelay(r) * 1e-3));
    
    % Packet drop rate (PDR)
%     dr1AvgPDR(r) = ( length(table2array(dr1Data(:, 6))) / (table2array(dr1Data(end, 6)) - table2array(dr1Data(1, 6))) ) * 100;
    dr2AvgPDR(r) = (1 - ( length(table2array(dr2Data(:, 6))) / ( table2array(dr2Data(end, 6)) - table2array(dr2Data(1, 6))+1 ) ) ) * 100;
    dr3AvgPDR(r) = (1 - ( length(table2array(dr3Data(:, 6))) / ( table2array(dr3Data(end, 6)) - table2array(dr3Data(1, 6))+1 ) ) ) * 100;
    dr4AvgPDR(r) = (1 - ( length(table2array(dr4Data(:, 6))) / ( table2array(dr4Data(end, 6)) - table2array(dr4Data(1, 6))+1 ) ) ) * 100;
    solarAvgPDR(r) = (1 - (length(table2array(solarData(:, 6))) / ( table2array(solarData(end, 6)) - table2array(solarData(1, 6))+1 )) ) * 100;
    pmuAvgPDR(r) = (1 - (length(table2array(pmuData(:, 6))) / ( table2array(pmuData(end, 6)) - table2array(pmuData(1, 6))+1 )) ) * 100;
end

% drDelay = mean([dr1AvgDelay; dr2AvgDelay; dr3AvgDelay; dr4AvgDelay]);
drDelay = mean([dr2AvgDelay; dr3AvgDelay; dr4AvgDelay]);
solarDelay = mean(solarAvgDelay);
pmuDelay = mean(pmuAvgDelay);

% drTpt = mean([dr1AvgTpt; dr2AvgTpt; dr3AvgTpt; dr4AvgTpt]);
drTpt = mean([dr2AvgTpt; dr3AvgTpt; dr4AvgTpt]);
solarTpt = mean(solarAvgTpt);
pmuTpt = mean(pmuAvgTpt);

% drPDR = mean([dr1AvgPDR; dr2AvgPDR; dr3AvgPDR; dr4AvgPDR]);
drPDR = mean([dr2AvgPDR; dr3AvgPDR; dr4AvgPDR]);
solarPDR = mean(solarAvgPDR);
pmuPDR = mean(pmuAvgPDR);

% Confidence intervals
drDelayCI = ConfInt([dr2AvgDelay; dr3AvgDelay; dr4AvgDelay]);
solarDelayCI = ConfInt(solarAvgDelay);
pmuDelayCI = ConfInt(pmuAvgDelay);

drPDRCI = ConfInt([dr2AvgPDR; dr3AvgPDR; dr4AvgPDR]);
solarPDRCI = ConfInt(solarAvgPDR);
pmuPDRCI = ConfInt(pmuAvgPDR);
% *********************************************************************** %

%% Convergence computations 
dr1AvgReward = [];
dr2AvgReward = [];
dr3AvgReward = [];
dr4AvgReward = [];
solarAvgReward = [];
pmuAvgReward = [];

dr2AvgCQI = [];
dr3AvgCQI = [];
dr4AvgCQI = [];
solarAvgCQI = [];
pmuAvgCQI = [];

for r = 1:nRuns
    convfileName = strcat('./BetaResults/convLogFile_', num2str(neNBs), 'eNBs_', num2str(nD2D),...
        'D2Ds_', num2str(nDummy), 'Dummy_', num2str(nUEs), 'UEs_', d2dScheduler, '_Run', num2str(r),...
        '_Beta', num2str(beta), '.dat');
    
    dataConv = importfileConv(convfileName);

    txID = dataConv{:, 2};
    
    dr1Data = dataConv(txID == 1, :);
    dr2Data = dataConv(txID == 2, :);
    dr3Data = dataConv(txID == 3, :);
    dr4Data = dataConv(txID == 4, :);
    solarData = dataConv(txID == 7, :);
    pmuData = dataConv(txID == 5, :);
        
    drLength = length(table2array(dr2Data(:, 4)));
    
    dr1Data = dr1Data((table2array(dr1Data(:, 3)) > explTime), :); 
    dr2Data = dr2Data((table2array(dr2Data(:, 3)) > explTime), :); 
    dr3Data = dr3Data((table2array(dr3Data(:, 3)) > explTime), :); 
    dr4Data = dr4Data((table2array(dr4Data(:, 3)) > explTime), :); 
    solarData = solarData((table2array(solarData(:, 3)) > explTime), :); 
    pmuData = pmuData((table2array(pmuData(:, 3)) > explTime), :); 
    
    % reward
%     dr1AvgReward(r) = table2array(dr1Data(:, 4));
    dr2AvgReward = [dr2AvgReward, table2array(dr2Data(:, 4))];
    dr3AvgReward = [dr3AvgReward, table2array(dr3Data(:, 4))];
    dr4AvgReward = [dr4AvgReward, table2array(dr4Data(:, 4))];
    solarAvgReward = [solarAvgReward, table2array(solarData(:, 4))];
    pmuAvgReward = [pmuAvgReward, table2array(pmuData(:, 4))];
    
    % CQI 
    dr2AvgCQI = [dr2AvgCQI, table2array(dr2Data(:, 5))];
    dr3AvgCQI = [dr3AvgCQI, table2array(dr3Data(:, 5))];
    dr4AvgCQI = [dr4AvgCQI, table2array(dr4Data(:, 5))];
    solarAvgCQI = [solarAvgCQI, table2array(solarData(:, 5))];
    pmuAvgCQI = [pmuAvgCQI, table2array(pmuData(:, 5))];
end

dr2Reward = mean(dr2AvgReward, 2);
dr3Reward = mean(dr3AvgReward, 2);
dr4Reward = mean(dr4AvgReward, 2);
drReward = mean([dr2Reward, dr3Reward, dr4Reward], 2);
solarReward = mean(solarAvgReward, 2);
pmuReward = mean(pmuAvgReward, 2);

drReward = cummean(drReward, 2);
solarReward = cummean(solarReward, 2);
pmuReward = cummean(pmuReward, 2);

% CQI average value
dr2CQI = mean(dr2AvgCQI, 2);
dr3CQI = mean(dr3AvgCQI, 2);
dr4CQI = mean(dr4AvgCQI, 2);
drCQI = mean([dr2CQI, dr3CQI, dr4CQI], 2);
solarCQI = mean(solarAvgCQI, 2);
pmuCQI = mean(pmuAvgCQI, 2);

drCQI = mean(drCQI);
solarCQI = mean(solarCQI);
pmuCQI = mean(pmuCQI);

% Confidence interval
drCQICI = ConfInt(mean([dr2CQI, dr3CQI, dr4CQI], 1));
solarCQICI = ConfInt(mean(solarAvgCQI, 1));
pmuCQICI = ConfInt(mean(pmuAvgCQI, 1));
% *********************************************************************** %


end








