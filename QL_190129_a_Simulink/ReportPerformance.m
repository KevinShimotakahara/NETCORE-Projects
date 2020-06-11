
function [SGD_AvgPktDelay, UE_AvgPktDelay, UNB_AvgPktDelay, SGD_Rate,...
    UE_Rate, UNB_Rate, SGD_maxRate, UE_maxRate, UNB_maxRate,...
    SGD_AvgWaitTime, UE_AvgWaitTime, UNB_AvgWaitTime,...
    SGD_RBs, UE_RBs, UNB_RBs, SGD_SFRBs, UE_SFRBs, UNB_SFRBs] = ...
    ReportPerformance(FileName, nSBSs, neNBs, nSGDs, nUEs, nUNBs)

global RunTime
DMQDelay = importfile(FileName);

%% SGDs
SGDsIdx             = find(DMQDelay{:, 1} == 'SGD');
MaxPktNumber        = max(DMQDelay{:, 5});
SGD_AvgPktDelay     = zeros(MaxPktNumber, nSGDs, nSBSs);
PktSize             = zeros(MaxPktNumber, nSGDs, nSBSs);
WaitTime            = zeros(MaxPktNumber, nSGDs, nSBSs);
SGD_RBs             = zeros(MaxPktNumber, nSGDs, nSBSs);
SGD_SFRBs           = zeros(RunTime, nSGDs, nSBSs);

SGD_DMQTable    = DMQDelay(SGDsIdx, :);
SGDs_Data       = table2array(SGD_DMQTable(:, 5));
SGDs_Tx         = table2array(SGD_DMQTable(:, 8));

for sbs_ = 1:nSBSs
    for s_ = 1:nSGDs
        for p_ = 1:MaxPktNumber
            SGDs_Data_Idx   = find(SGDs_Data == p_);
            SGDs_Data1      = SGD_DMQTable(SGDs_Data_Idx, :);
            SGDs_Data2      = SGDs_Data1((SGDs_Data1{:,2} == s_), :);
            SGDs_Data3      = SGDs_Data2((SGDs_Data2{:, 4} == sbs_), :);
            PktsDelays      = SGDs_Data3{:, 9} - SGDs_Data3{:, 7};
            SGD_AvgPktDelay(p_, s_, sbs_)   = sum(PktsDelays);
            PktSize_temp                    = SGDs_Data3{:, 6};
            PktSize(p_, s_, sbs_)           = sum(PktSize_temp);
            WaitTime(p_, s_, sbs_)          = sum(SGDs_Data3{:, 8} - SGDs_Data3{:, 7});
            SGD_RBs(p_, s_, sbs_)           = sum(SGDs_Data3{:, 10});
        end
        
        for sf_ = 1:RunTime
            SGDs_Data_Idx   = find(SGDs_Tx == sf_);
            SGDs_Data1      = SGD_DMQTable(SGDs_Data_Idx, :);
            SGDs_Data2      = SGDs_Data1((SGDs_Data1{:,2} == s_), :);
            SGDs_Data3      = SGDs_Data2((SGDs_Data2{:, 4} == sbs_), :);
            if(~isempty(SGDs_Data3))
                SGD_SFRBs(sf_, s_, sbs_)  = SGDs_Data3{1, 10};
            end
        end
    end
end

SGD_AvgPktDelay_temp = [];
PktSize_temp         = [];
WaitTime_temp        = [];
SGD_RBs_temp         = [];
SGD_SFRBs_temp       = [];

for j = 1:size(SGD_AvgPktDelay, 3)
    for i = 1:size(SGD_AvgPktDelay, 2)
        SGD_AvgPktDelay_temp = [SGD_AvgPktDelay_temp, (SGD_AvgPktDelay(:, i, j))];
        PktSize_temp         = [PktSize_temp, (PktSize(:, i, j))];
        WaitTime_temp        = [WaitTime_temp, WaitTime(:, i, j)];
        SGD_RBs_temp         = [SGD_RBs_temp, SGD_RBs(:, i, j)];
        SGD_SFRBs_temp       = [SGD_SFRBs_temp, SGD_SFRBs(:, i, j)];
    end
end

% Delay
SGD_AvgPktDelay = SGD_AvgPktDelay_temp;
PktSize         = PktSize_temp;

% Throughput
SGD_Rate    = PktSize ./ (SGD_AvgPktDelay * 1e-3);
SGD_Rate(isnan(SGD_Rate)) = 0;
SGD_Rate    = sum(SGD_Rate(end-15:end, :), 1);
SGD_maxRate = max(SGD_Rate) / 1e6;
SGD_Rate    = mean(SGD_Rate);

% Average waiting time
SGD_AvgWaitTime         = WaitTime_temp;

% RBs to SGDs
SGD_RBs     = SGD_RBs_temp;
SGD_SFRBs   = SGD_SFRBs_temp; 
% *********************************************************************** %

%% UEs
SGDsIdx             = find(DMQDelay{:, 1} == 'UE');
MaxPktNumber        = max(DMQDelay{:, 5});
UE_AvgPktDelay      = zeros(MaxPktNumber, nUEs, nSBSs);
PktSize             = zeros(MaxPktNumber, nUEs, nSBSs);
WaitTime            = zeros(MaxPktNumber, nUEs, nSBSs);
UE_RBs              = zeros(MaxPktNumber, nUEs, nSBSs);
UE_SFRBs            = zeros(RunTime, nSGDs, nSBSs);

SGD_DMQTable    = DMQDelay(SGDsIdx, :);
SGDs_Data       = table2array(SGD_DMQTable(:, 5));
SGDs_Tx         = table2array(SGD_DMQTable(:, 8));
            
for sbs_ = 1:nSBSs
    for s_ = 1:nUEs
        for p_ = 1:MaxPktNumber
            UE_id           = s_ + nSGDs;
            SGDs_Data_Idx   = find(SGDs_Data == p_);
            SGDs_Data1      = SGD_DMQTable(SGDs_Data_Idx, :);
            SGDs_Data2     = SGDs_Data1((SGDs_Data1{:, 2} == UE_id), :);
            SGDs_Data3     = SGDs_Data2((SGDs_Data2{:, 4} == sbs_), :);
            PktsDelays     = SGDs_Data3{:, 9} - SGDs_Data3{:, 7};
            UE_AvgPktDelay(p_, s_, sbs_) = sum(PktsDelays);
            PktSize_temp                 = SGDs_Data3{:, 6};
            PktSize(p_, s_, sbs_)        = sum(PktSize_temp);
            WaitTime(p_, s_, sbs_)       = sum(SGDs_Data3{:, 8} - SGDs_Data3{:, 7});
            UE_RBs(p_, s_, sbs_)         = sum(SGDs_Data3{:, 10});
        end
        
        for sf_ = 1:RunTime
            UE_id           = s_ + nSGDs;
            SGDs_Data_Idx   = find(SGDs_Tx == sf_);
            SGDs_Data1      = SGD_DMQTable(SGDs_Data_Idx, :);
            SGDs_Data2      = SGDs_Data1((SGDs_Data1{:,2} == UE_id), :);
            SGDs_Data3      = SGDs_Data2((SGDs_Data2{:, 4} == sbs_), :);
            if(~isempty(SGDs_Data3))
                UE_SFRBs(sf_, s_, sbs_)  = SGDs_Data3{1, 10};
            end
        end
    end
end

UE_AvgPktDelay_temp = [];
PktSize_temp        = [];
WaitTime_temp       = [];
UE_RBs_temp         = [];
UE_SFRBs_temp       = [];

for j = 1:size(UE_AvgPktDelay, 3)
    for i = 1:size(UE_AvgPktDelay, 2)
        UE_AvgPktDelay_temp  = [UE_AvgPktDelay_temp, (UE_AvgPktDelay(:, i, j))];
        PktSize_temp         = [PktSize_temp, (PktSize(:, i, j))];
        WaitTime_temp        = [WaitTime_temp, WaitTime(:, i, j)];
        UE_RBs_temp          = [UE_RBs_temp, UE_RBs(:, i, j)];
        UE_SFRBs_temp        = [UE_SFRBs_temp, UE_SFRBs(:, i, j)];
    end
end

% Delay
UE_AvgPktDelay = UE_AvgPktDelay_temp;
PktSize         = PktSize_temp;

% Throughput
UE_Rate    = PktSize ./ (UE_AvgPktDelay * 1e-3);
UE_Rate(isnan(UE_Rate)) = 0;
UE_Rate    = sum(UE_Rate(end-15:end, :), 1);
UE_maxRate = max(UE_Rate) / 1e6;
UE_Rate    = mean(UE_Rate);

% Average waiting time
UE_AvgWaitTime            = WaitTime_temp;

% UEs RBs
UE_RBs      = UE_RBs_temp;
UE_SFRBs    = UE_SFRBs_temp;

%% UNBs
UNBsIdx             = find(DMQDelay{:, 1} == 'UNB');
MaxPktNumber        = max(DMQDelay{:, 5});
UNB_AvgPktDelay     = zeros(MaxPktNumber, nUNBs, neNBs);
PktSize             = zeros(MaxPktNumber, nUNBs, neNBs);
WaitTime            = zeros(MaxPktNumber, nUNBs, neNBs);
UNB_RBs             = zeros(MaxPktNumber, nUNBs, neNBs);
UNB_SFRBs           = zeros(RunTime, nSGDs, nSBSs);

SGD_DMQTable    = DMQDelay(UNBsIdx, :);
SGDs_Data       = table2array(SGD_DMQTable(:, 5));
SGDs_Tx         = table2array(SGD_DMQTable(:, 8));

for e_ = 1:neNBs
    for s_ = 1:nUNBs
        for p_ = 1:MaxPktNumber
            UNB_id          = s_;
            SGDs_Data_Idx   = find(SGDs_Data == p_);
            SGDs_Data1      = SGD_DMQTable(SGDs_Data_Idx, :);
            SGDs_Data2      = SGDs_Data1((SGDs_Data1{:, 2} == UNB_id), :);
            SGDs_Data3      = SGDs_Data2((SGDs_Data2{:, 4} == e_), :);
            PktsDelays      = SGDs_Data3{:, 9} - SGDs_Data3{:, 7};
            UNB_AvgPktDelay(p_, s_, e_) = sum(PktsDelays);
            PktSize_temp                 = SGDs_Data3{:, 6};
            PktSize(p_, s_, e_)        = sum(PktSize_temp);
            WaitTime(p_, s_, e_)       = sum(SGDs_Data3{:, 8} - SGDs_Data3{:, 7});
            UNB_RBs(p_, s_, e_)        = sum(SGDs_Data3{:, 10});
        end
        
        for sf_ = 1:RunTime
            UNB_id          = s_;
            SGDs_Data_Idx   = find(SGDs_Tx == sf_);
            SGDs_Data1      = SGD_DMQTable(SGDs_Data_Idx, :);
            SGDs_Data2      = SGDs_Data1((SGDs_Data1{:,2} == UNB_id), :);
            SGDs_Data3      = SGDs_Data2((SGDs_Data2{:, 4} == e_), :);
            if(~isempty(SGDs_Data3))
                UNB_SFRBs(sf_, s_, e_)  = SGDs_Data3{1, 10};
            end
        end
    end
end

UNB_AvgPktDelay_temp = [];
PktSize_temp         = [];
WaitTime_temp        = [];
UNB_RBs_temp         = [];
UNB_SFRBs_temp       = [];

for j = 1:size(UNB_AvgPktDelay, 3)
    for i = 1:size(UNB_AvgPktDelay, 2)
        UNB_AvgPktDelay_temp  = [UNB_AvgPktDelay_temp, (UNB_AvgPktDelay(:, i, j))];
        PktSize_temp         = [PktSize_temp, (PktSize(:, i, j))];
        WaitTime_temp        = [WaitTime_temp, WaitTime(:, i, j)];
        UNB_RBs_temp         = [UNB_RBs_temp, UNB_RBs(:, i, j)];
        UNB_SFRBs_temp       = [UNB_SFRBs_temp, UNB_SFRBs(:, i, j)];
    end
end

% Delay
UNB_AvgPktDelay = UNB_AvgPktDelay_temp;
PktSize         = PktSize_temp;

% Throughput
UNB_Rate    = PktSize ./ (UNB_AvgPktDelay * 1e-3);
UNB_Rate(isnan(UNB_Rate)) = 0;
UNB_Rate    = sum(UNB_Rate(end-15:end, :), 1);
UNB_maxRate = max(UNB_Rate) / 1e6;
UNB_Rate    = mean(UNB_Rate);

% Average waiting time
UNB_AvgWaitTime            = WaitTime_temp;

% UEs RBs
UNB_RBs     = UNB_RBs_temp;
UNB_SFRBs   = UNB_SFRBs_temp;
end



