
clear
clc

Seed = 2^3;
rng(Seed);

simTime = 1000;
nAgents = 6;

global gClk
gClk = 0;
% *********************************************************************** %

%% Log information
global logData
global convLogFile

logData = logdataClass;
convLogFile = './Results/logFile2.dat';

logData.report_fn(convLogFile);

convInfo.nodeID = 'nodeID,';
convInfo.nSF = 'nSF,';
convInfo.reward = 'reward,';
convInfo.cqi = 'cqi,';
convInfo.remDelay = 'remDelay,';
convInfo.prState = 'prState,';
convInfo.newState = 'newState,';
convInfo.action = 'action,';
convInfo.qValue = 'qValue,';
convInfo.epsilon = 'epsilon\r\n';

logData.saveconvLog_fn(convInfo);
% *********************************************************************** %

for a = 1:nAgents
    Agents(a) = agent;
    
    Agents(a).nodeID = a;
    Agents(a).oldest_crtTS = 0;
    Agents(a).cqiBuffer = 9;
    
    Agents(a).schedInfo.prAction = 1;
    Agents(a).schedInfo.prState = 1;
    
    Agents(a).scheduler = Schedulers.QLearning;
    Agents(a).scheduler.initScheduler_fn(10, 6, 30, 50, 8, 9, Agents(a), 0.5);
end

actIdx = zeros(1, nAgents);

while(gClk < simTime)
    for a = 1:nAgents
        [actIdx(a)] = Agents(a).scheduler.schedule_fn();
    end

    % punished agents
    for i = 1:6
       punishAgents = find(actIdx == i);
       if(length(punishAgents) >= 2)
           for j = punishAgents
               Agents(j).oldest_crtTS = 60;
               Agents(j).cqiBuffer = 1;
           end
       end
    end
    
    % rewarded agents
    actIdxOccur = tabulate(actIdx);
    rewardedActions = actIdxOccur(actIdxOccur(:, 2) == 1)';
    for i = rewardedActions
        rewardAgent = find(actIdx == i);
        Agents(rewardAgent).oldest_crtTS = 0;
        Agents(rewardAgent).cqiBuffer = 9;
    end
           
    gClk = gClk + 1;
    
    fprintf('Simulation time = %d\n', gClk);
end

























