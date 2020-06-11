


nRuns = 1;
Dummy = 0;
beta = 0.5;
explTime = 18000;        % time after which results are computed

nDummy = length(Dummy);

%% Q-Learning results
d2dScheduler = 'QLearning';

i = 1;

QLdrDelay = zeros(nDummy, 1);
QLsolarDelay = zeros(nDummy, 1);
QLpmuDelay = zeros(nDummy, 1);

QLdrTpt = zeros(nDummy, 1);
QLsolarTpt = zeros(nDummy, 1);
QLpmuTpt = zeros(nDummy, 1);

QLdrPDR = zeros(nDummy, 1);
QLsolarPDR = zeros(nDummy, 1);
QLpmuPDR = zeros(nDummy, 1);

QLdrCQI = zeros(nDummy, 1);
QLsolarCQI = zeros(nDummy, 1);
QLpmuCQI = zeros(nDummy, 1);

QLdrDelayCI = [];
QLsolarDelayCI = [];
QLpmuDelayCI = [];
QLdrPDRCI = []; 
QLsolarPDRCI = []; 
QLpmuPDRCI = [];
QLdrCQICI = [];
QLsolarCQICI = [];
QLpmuCQICI = [];
    
for b = Dummy
    [QLdrDelay(i), QLsolarDelay(i), QLpmuDelay(i),...
        QLdrTpt(i), QLsolarTpt(i), QLpmuTpt(i),...
        QLdrPDR(i), QLsolarPDR(i), QLpmuPDR(i),...   
        QLdrReward, QLsolarReward, QLpmuReward,...
        QLdrCQI(i), QLsolarCQI(i), QLpmuCQI(i),...
        QLdrDelayCITemp, QLsolarDelayCITemp, QLpmuDelayCITemp,...
        QLdrPDRCITemp, QLsolarPDRCITemp, QLpmuPDRCITemp,...
        QLdrCQICITemp, QLsolarCQICITemp, QLpmuCQICITemp] =...
    extractResults(nRuns, b, beta, d2dScheduler, explTime);  
    
    QLdrDelayCI = [QLdrDelayCI; QLdrDelayCITemp];
    QLsolarDelayCI = [QLsolarDelayCI; QLsolarDelayCITemp];
    QLpmuDelayCI = [QLpmuDelayCI; QLpmuDelayCITemp];
    
    QLdrPDRCI = [QLdrPDRCI; QLdrPDRCITemp]; 
    QLsolarPDRCI = [QLsolarPDRCI; QLsolarPDRCITemp]; 
    QLpmuPDRCI = [QLpmuPDRCI; QLpmuPDRCITemp];
    
    QLdrCQICI = [QLdrCQICI; QLdrCQICITemp];
    QLsolarCQICI = [QLsolarCQICI; QLsolarCQICITemp];
    QLpmuCQICI = [QLpmuCQICI; QLpmuCQICITemp];
    
    i = i+1;
end
% *********************************************************************** %

