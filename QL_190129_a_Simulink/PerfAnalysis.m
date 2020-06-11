
function [SGD_TotAvgDelay, UE_TotAvgDelay, UNB_TotAvgDelay,...
        SGD_TotAvgRate, UE_TotAvgRate, UNB_TotAvgRate,...
        SGD_Outage, UE_Outage, UNB_Outage,...
        JFI, SGD_TopPRate, UE_TopPRate, UNB_TopPRate,...
        SGD_AvgWtTime, UE_AvgWtTime, UNB_AvgWtTime, ...
        SGD_RBs, UE_RBs, UNB_RBs] = ...
    PerfAnalysis(SGD_PktDelay, UE_PktDelay, UNB_PktDelay,... 
        SGD_Rate, UE_Rate, UNB_Rate,...
        SGD_AvgWaitTime, UE_AvgWaitTime, UNB_AvgWaitTime, ...
        SGD_RBs, UE_RBs, UNB_RBs)

global nSGDs
global nUEs
global nSBSs

StartCapture    = 1; % 15
EndCapture      = find(min(SGD_PktDelay, [], 2) == 0, 1);   % 20
EvalPer         = 20; 

%% Average delay and rate
% Packet delay
SGD_AvgPktDelay_new = SGD_PktDelay(end-EvalPer:end, :);
SGD_AvgPktDelay_new = mean(SGD_AvgPktDelay_new, 2);
SGD_TotAvgDelay     = mean(SGD_AvgPktDelay_new);

UE_AvgPktDelay_new  = UE_PktDelay(end-EvalPer:end, :);
UE_AvgPktDelay_new  = mean(UE_AvgPktDelay_new, 2);
UE_TotAvgDelay      = mean(UE_AvgPktDelay_new);

UNB_AvgPktDelay_new  = UNB_PktDelay(end-EvalPer:end, :);
UNB_AvgPktDelay_new  = mean(UNB_AvgPktDelay_new, 2);
UNB_TotAvgDelay      = mean(UNB_AvgPktDelay_new);

% Wait time
SGD_WtTime_new      = SGD_AvgWaitTime(end-EvalPer:end, :);
SGD_WtTime_new      = mean(SGD_WtTime_new, 2);
SGD_AvgWtTime       = mean(SGD_WtTime_new);

UE_WtTime_new       = UE_AvgWaitTime(end-EvalPer:end, :);
UE_WtTime_new       = mean(UE_WtTime_new, 2);
UE_AvgWtTime        = mean(UE_WtTime_new);

UNB_WtTime_new       = UNB_AvgWaitTime(end-EvalPer:end, :);
UNB_WtTime_new       = mean(UNB_WtTime_new, 2);
UNB_AvgWtTime        = mean(UNB_WtTime_new);

% Rate
SGD_TotAvgRate = sum(SGD_Rate) / 1e6;
UE_TotAvgRate = sum(UE_Rate) / 1e6;
UNB_TotAvgRate = sum(UNB_Rate) / 1e6;
% ----------------------------------------------------------------------- %

%% Outage
All_Delay = [SGD_PktDelay];
All_Delay = mean(All_Delay(StartCapture:end-EndCapture, :), 1);
Outage = All_Delay(All_Delay == 0);
Outage = length(Outage);
SGD_Outage = (Outage / (nSBSs*nSGDs))*100;

All_Delay = [UE_PktDelay];
All_Delay = mean(All_Delay(StartCapture:end-EndCapture, :), 1);
Outage = All_Delay(All_Delay == 0);
Outage = length(Outage);
UE_Outage = (Outage / (nSBSs*nUEs))*100;

All_Delay = [UNB_PktDelay];
All_Delay = mean(All_Delay(StartCapture:end-EndCapture, :), 1);
Outage = All_Delay(All_Delay == 0);
Outage = length(Outage);
UNB_Outage = (Outage / (nSBSs*nUEs))*100;
% *********************************************************************** %

%% Spectral efficiency
SGD_RBs     = sum(SGD_RBs, 1);
UE_RBs      = sum(UE_RBs, 1);
UNB_RBs     = sum(UNB_RBs, 1);

% *********************************************************************** %

%% Jain's fairness index
nRBs = [SGD_RBs, UE_RBs, UNB_RBs];
JFI = (sum(nRBs))^2 / (length(nRBs) * sum(nRBs.^2));
% *********************************************************************** %

%% Entropy fairness
nRBs(nRBs == 0) = [];
Psum = sum(nRBs);
P = nRBs / Psum;

EFI = sum(P .* log2(1./P));
% *********************************************************************** %

%% Top-10 max throughput users
SGD_TopPRate = sort(SGD_Rate, 'descend');

if(length(SGD_TopPRate) > 10)
    SGD_TopPRate = mean(SGD_TopPRate(1:10));
else
    SGD_TopPRate = mean(SGD_TopPRate);
end

SGD_TopPRate = SGD_TopPRate/1e6;

UE_TopPRate = sort(UE_Rate, 'descend');

if(length(UE_TopPRate) > 10)
    UE_TopPRate = mean(UE_TopPRate(1:10));
else
    UE_TopPRate = mean(UE_TopPRate);
end

UE_TopPRate = UE_TopPRate/1e6;

UNB_TopPRate = sort(UNB_Rate, 'descend');

if(length(UNB_TopPRate) > 10)
    UNB_TopPRate = mean(UNB_TopPRate(1:10));
else
    UNB_TopPRate = mean(UNB_TopPRate);
end

UNB_TopPRate = UNB_TopPRate/1e6;

% *********************************************************************** %


end





