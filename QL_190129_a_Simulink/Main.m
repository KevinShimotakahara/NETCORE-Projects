
%% Main function: 
% This is the main function of the simulator, it is used for ul/sl
% transmission/reception. 

function Main()

%% Global variables
global gClk    
global nRun
global RunTime
global Seed

global eNBs             % structure for eNB
global UDs              % structure for user device

%% Display subframe number
fprintf('\n******************************************************\n');
fprintf('**************** [Subframe number = %d] ***************\n', gClk);
fprintf('******************************************************\n');
% *************************************************************** %

%% Simulation main loop
% Steps of simulation (TTI-based)
% 1. UEs: Check txBuffer for data to transmit.
% 2. Perform scheduling on D2D devices.
% 3. UDs send SR or perform sl/ul tx. D2D perform sidelink scheduling.
% 4. eNB rx ul tx or SR from the users.
% 5. eNB tx (on dl) pools configuration to D2D or ul grants to UEs.
% 6. Check the SR and update scheduler buffers.
% *************************************************************** %

%% 1. UDs/UEs: Check txBuffer for data to transmit.
% D2D get their traffic from the power circuit (No action required here).
% UEs generate their traffic using the poisson distribution.
for u = 1:length(UDs)
    UDs(u).trafficGen_fn();
end
% *************************************************************** %

%% 3. UDs send SR or transmit
% UEs transmit SR whenever there is new data, or ul transmit when there is
% ul grant.
% D2Ds perform sl tx.
% actIdx = [];
for u = 1:length(UDs)
    UDs(u).tx_fn();
%     actIdx = [actIdx, UDs(u).schedInfo.prAction];
end

% punished agents
% for i = 1:7
%    punishAgents = find(actIdx == i);
%    if(length(punishAgents) >= 2)
%        for j = punishAgents
%            UDs(j).cqiBuffer = 1;
%        end
%    end
% end

% rewarded agents
% actIdxOccur = tabulate(actIdx);
% rewardedActions = actIdxOccur(actIdxOccur(:, 2) == 1)';
% for i = rewardedActions
%     UDs(actIdx == i).cqiBuffer = 9;
% end
% *************************************************************** %

%% 4. eNB rx ul tx or SR from the users
for e = 1:length(eNBs)
    eNBs(e).ulrx_fn();
end
% *************************************************************** %

%% 5. D2D rx on sidelink
for u = 1:length(UDs)
    UDs(u).slrx_fn();
end
% *************************************************************** %

%% 6. eNB tx (on dl) pools configuration to D2D or ul grants to UEs
% in both cases, this is performed 4 subframes after receiving the requests
% from the users.
% Note that: D2D pool configuration is updated every PSCCH period only.
% This is triggered by D2D SR request at the end of the PSCCH period.
for e = 1:length(eNBs)
    eNBs.dltx_fn();
end
% *************************************************************** %

%% 7. Check the SR and update scheduler buffers. And check for pools configuration
for e = 1:length(eNBs)
    eNBs(e).checkSR_fn();
end
% *************************************************************** %

%% Advance the clock by 1 (1 TTI)
gClk = gClk + 1;
% *************************************************************** %

%% Reset for each Run
if(gClk >= RunTime)
    nRun = nRun + 1;
    gClk = 0;
    Seed = 2*Seed;
    rng(Seed);
    Main_Init();
end
% *************************************************************** %

end







