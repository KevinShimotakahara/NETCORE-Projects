
clear
clc

nRuns = 15;
Dummy = [0 1 2 3 4 5];
beta = 0;
explTime = 10000;        % time after which results are computed

nDummy = length(Dummy);

%% Random results
d2dScheduler = 'RandomNew';

i = 1;

RandomdrDelay = zeros(nDummy, 1);
RandomsolarDelay = zeros(nDummy, 1);
RandompmuDelay = zeros(nDummy, 1);

RandomdrTpt = zeros(nDummy, 1);
RandomsolarTpt = zeros(nDummy, 1);
RandompmuTpt = zeros(nDummy, 1);

RandomdrPDR = zeros(nDummy, 1);
RandomsolarPDR = zeros(nDummy, 1);
RandompmuPDR = zeros(nDummy, 1);

RandomdrCQI = zeros(nDummy, 1);
RandomsolarCQI = zeros(nDummy, 1);
RandompmuCQI = zeros(nDummy, 1);

RandomdrDelayCI = [];
RandomsolarDelayCI = [];
RandompmuDelayCI = [];
RandomdrPDRCI = []; 
RandomsolarPDRCI = []; 
RandompmuPDRCI = [];
RandomdrCQICI = [];
RandomsolarCQICI = [];
RandompmuCQICI = [];
    
for b = Dummy
    [RandomdrDelay(i), RandomsolarDelay(i), RandompmuDelay(i),...
        RandomdrTpt(i), RandomsolarTpt(i), RandompmuTpt(i),...
        RandomdrPDR(i), RandomsolarPDR(i), RandompmuPDR(i),...   
        RandomdrReward, RandomsolarReward, RandompmuReward,...
        RandomdrCQI(i), RandomsolarCQI(i), RandompmuCQI(i),...
        RandomdrDelayCITemp, RandomsolarDelayCITemp, RandompmuDelayCITemp,...
        RandomdrPDRCITemp, RandomsolarPDRCITemp, RandompmuPDRCITemp,...
        RandomdrCQICITemp, RandomsolarCQICITemp, RandompmuCQICITemp] =...
    extractResults(nRuns, b, beta, d2dScheduler, explTime);  
    
    RandomdrDelayCI = [RandomdrDelayCI; RandomdrDelayCITemp];
    RandomsolarDelayCI = [RandomsolarDelayCI; RandomsolarDelayCITemp];
    RandompmuDelayCI = [RandompmuDelayCI; RandompmuDelayCITemp];
    
    RandomdrPDRCI = [RandomdrPDRCI; RandomdrPDRCITemp]; 
    RandomsolarPDRCI = [RandomsolarPDRCI; RandomsolarPDRCITemp]; 
    RandompmuPDRCI = [RandompmuPDRCI; RandompmuPDRCITemp];
    
    RandomdrCQICI = [RandomdrCQICI; RandomdrCQICITemp];
    RandomsolarCQICI = [RandomsolarCQICI; RandomsolarCQICITemp];
    RandompmuCQICI = [RandompmuCQICI; RandompmuCQICITemp];
    
    i = i+1;
end
% *********************************************************************** %

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

%% Average delay
figure; 
subplot(3, 1, 1);
grid on; hold on
h1 = errorbar(Dummy, RandomdrDelay, (RandomdrDelay-RandomdrDelayCI(:, 1)), (RandomdrDelayCI(:, 2)-RandomdrDelay), 'r');
h1 = [h1; errorbar(Dummy, QLdrDelay, (QLdrDelay-QLdrDelayCI(:, 1)), (QLdrDelayCI(:, 2)-QLdrDelay), 'b')];
xlabel('Number of dummy devices');
ylabel('DR delay [msec]');
legend('Random', 'Q-Learning');

subplot(3, 1, 2);
grid on; hold on
h2 = errorbar(Dummy, RandomsolarDelay, (RandomsolarDelay-RandomsolarDelayCI(:, 1)), (RandomsolarDelayCI(:, 2)-RandomsolarDelay), 'r');
h2 = [h2; errorbar(Dummy, QLsolarDelay, (QLsolarDelay-QLsolarDelayCI(:, 1)), (QLsolarDelayCI(:, 2)-QLsolarDelay), 'b')];
xlabel('Number of dummy devices');
ylabel('Solar delay [msec]');
legend('Random', 'Q-Learning');

subplot(3, 1, 3);
grid on; hold on
h3 = errorbar(Dummy, RandompmuDelay, (RandompmuDelay-RandompmuDelayCI(:, 1)), (RandompmuDelayCI(:, 2)-RandompmuDelay), 'r');
h3 = [h3; errorbar(Dummy, QLpmuDelay, (QLpmuDelay-QLpmuDelayCI(:, 1)), (QLpmuDelayCI(:, 2)-QLpmuDelay), 'b')];
xlabel('Number of dummy devices');
ylabel('PMU delay [msec]');
legend('Random', 'Q-Learning');

gca.FontSize = 14;
gca.FontWeight = 'Bold';
set(gcf, 'Position', [400, 200, 1000, 700]);
% *********************************************************************** %

%% PDR
figure; 
subplot(3, 1, 1);
grid on; hold on
h1 = errorbar(Dummy, RandomdrPDR, (RandomdrPDR-RandomdrPDRCI(:, 1)), (RandomdrPDRCI(:, 2)-RandomdrPDR), 'r');
h1 = [h1; errorbar(Dummy, QLdrPDR, (QLdrPDR-QLdrPDRCI(:, 1)), (QLdrPDRCI(:, 2)-QLdrPDR), 'b')];
xlabel('Number of dummy devices');
ylabel('DR PDR [msec]');
legend('Random', 'Q-Learning');

subplot(3, 1, 2);
grid on; hold on
h2 = errorbar(Dummy, RandomsolarPDR, (RandomsolarPDR-RandomsolarPDRCI(:, 1)), (RandomsolarPDRCI(:, 2)-RandomsolarPDR), 'r');
h2 = [h2; errorbar(Dummy, QLsolarPDR, (QLsolarPDR-QLsolarPDRCI(:, 1)), (QLsolarPDRCI(:, 2)-QLsolarPDR), 'b')];
xlabel('Number of dummy devices');
ylabel('Solar PDR [%]');
legend('Random', 'Q-Learning');

subplot(3, 1, 3);
grid on; hold on
h3 = errorbar(Dummy, RandompmuPDR, (RandompmuPDR-RandompmuPDRCI(:, 1)), (RandompmuPDRCI(:, 2)-RandompmuPDR), 'r');
h3 = [h3; errorbar(Dummy, QLpmuPDR, (QLpmuPDR-QLpmuPDRCI(:, 1)), (QLpmuPDRCI(:, 2)-QLpmuPDR), 'b')];
xlabel('Number of dummy devices');
ylabel('PMU PDR [%]');
legend('Random', 'Q-Learning');

gca.FontSize = 14;
gca.FontWeight = 'Bold';
set(gcf, 'Position', [400, 200, 1000, 700]);
% *********************************************************************** %

%% CQI
% figure; 
% subplot(3, 1, 1);
% grid on; hold on
% h1 = errorbar(Dummy, RandomdrCQI, (RandomdrPDR-RandomdrCQICI(:, 1)), (RandomdrCQICI(:, 2)-RandomdrCQI), 'r');
% h1 = [h1; errorbar(Dummy, QLdrCQI, (QLdrCQI-QLdrCQICI(:, 1)), (QLdrCQICI(:, 2)-QLdrCQI), 'b')];
% xlabel('Number of dummy devices');
% ylabel('DR CQI [msec]');
% legend('Random', 'Q-Learning');
% 
% subplot(3, 1, 2);
% grid on; hold on
% h2 = errorbar(Dummy, RandomsolarCQI, (RandomsolarCQI-RandomsolarCQICI(:, 1)), (RandomsolarCQICI(:, 2)-RandomsolarCQI), 'r');
% h2 = [h2; errorbar(Dummy, QLsolarCQI, (QLsolarCQI-QLsolarCQICI(:, 1)), (QLsolarCQICI(:, 2)-QLsolarCQI), 'b')];
% xlabel('Number of dummy devices');
% ylabel('Solar CQI [msec]');
% legend('Random', 'Q-Learning');
% 
% subplot(3, 1, 3);
% grid on; hold on
% h3 = errorbar(Dummy, RandompmuCQI, (RandompmuCQIy-RandompmuCQICI(:, 1)), (RandompmuCQICI(:, 2)-RandompmuCQI), 'r');
% h3 = [h3; errorbar(Dummy, QLpmuCQI, (QLpmuDelay-QLpmuCQICI(:, 1)), (QLpmuCQICI(:, 2)-QLpmuCQI), 'b')];
% xlabel('Number of dummy devices');
% ylabel('PMU CQI [msec]');
% legend('Random', 'Q-Learning');
% 
% gca.FontSize = 14;
% gca.FontWeight = 'Bold';
% set(gcf, 'Position', [400, 200, 1000, 700]);
% *********************************************************************** %


