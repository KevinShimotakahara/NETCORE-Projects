
clear
clc

%% Plot results
global d2dScheduler
d2dScheduler = 'QLearning';

global neNBs            % number of eNBs
global nD2D             % number of SGDs
global nDummy
global nUEs             % number of UEs

neNBs = 1;          % Number of eNBs
nD2D = 7;           % Number of D2D users per eNB
nDummy = 4;
nUEs = 0;           % Number of UEs users per eNB

nRuns = 10;
beta = 0:0.25:0.75;


i = 1;
nbeta = length(beta);
drDelay = zeros(nbeta, 1);
solarDelay = zeros(nbeta, 1);
pmuDelay = zeros(nbeta, 1);

drTpt = zeros(nbeta, 1);
solarTpt = zeros(nbeta, 1);
pmuTpt = zeros(nbeta, 1);

drPDR = zeros(nbeta, 1);
solarPDR = zeros(nbeta, 1);
pmuPDR = zeros(nbeta, 1);

drCQI = zeros(nbeta, 1);
solarCQI = zeros(nbeta, 1);
pmuCQI = zeros(nbeta, 1);

drDelayCI = [];
solarDelayCI = [];
pmuDelayCI = [];
drPDRCI = []; 
solarPDRCI = []; 
pmuPDRCI = [];
drCQICI = [];
solarCQICI = [];
pmuCQICI = [];

for b = beta
    [drDelay(i), solarDelay(i), pmuDelay(i),...
        drTpt(i), solarTpt(i), pmuTpt(i),...
        drPDR(i), solarPDR(i), pmuPDR(i),...   
        drReward, solarReward, pmuReward,...
        drCQI(i), solarCQI(i), pmuCQI(i),...
        drDelayCITemp, solarDelayCITemp, pmuDelayCITemp,...
        drPDRCITemp, solarPDRCITemp, pmuPDRCITemp,...
        drCQICITemp, solarCQICITemp, pmuCQICITemp] =...
    extractBetaResults(nRuns, b);
        
    drDelayCI = [drDelayCI; drDelayCITemp];
    solarDelayCI = [solarDelayCI; solarDelayCITemp];
    pmuDelayCI = [pmuDelayCI; pmuDelayCITemp];
    
    drPDRCI = [drPDRCI; drPDRCITemp]; 
    solarPDRCI = [solarPDRCI; solarPDRCITemp]; 
    pmuPDRCI = [pmuPDRCI; pmuPDRCITemp];
    
    drCQICI = [drCQICI; drCQICITemp];
    solarCQICI = [solarCQICI; solarCQICITemp];
    pmuCQICI = [pmuCQICI; pmuCQICITemp];  
    
    i = i+1;
end
% *********************************************************************** %

%% Plot delay
% figure; 
% h1 = errorbar(beta, drDelay, (drDelay-drDelayCI(:, 1)), (drDelayCI(:, 2)-drDelay), 'b');
% grid on; hold on
% h1 = [h1; errorbar(beta, solarDelay, (solarDelay-solarDelayCI(:, 1)), (solarDelayCI(:, 2)-solarDelay), 'r')];
% h1 = [h1; errorbar(beta, pmuDelay, (pmuDelay-pmuDelayCI(:, 1)), (pmuDelayCI(:, 2)-pmuDelay), 'k')];
% xlabel('PMU \beta');
% ylabel('Average packet delay [msec]');
% legend('DR', 'Solar', 'PMU', 'location', 'best');
% ax = gca;
% ax.FontSize = 14;
% ax.FontWeight = 'Bold';
% set(h1, 'LineWidth', 2);
% % xlim([0 1]);
% % ylim([6 25]);
% set(gcf, 'Position', [500, 300, 800, 500]);
% *********************************************************************** %

%% Plot PDR
% figure; 
% h1 = errorbar(beta, drPDR, (drPDR-drPDRCI(:, 1)), (drPDRCI(:, 2)-drPDR), 'b');
% grid on; hold on
% h1 = [h1; errorbar(beta, solarPDR, (solarPDR-solarPDRCI(:, 1)), (solarPDRCI(:, 2)-solarPDR), 'r')];
% h1 = [h1; errorbar(beta, pmuPDR, (pmuPDR-pmuPDRCI(:, 1)), (pmuPDRCI(:, 2)-pmuPDR), 'k')];
% % title('Average packet delay');
% xlabel('PMU \beta');
% ylabel('Packet Drop Rate (PDR)');
% legend('DR', 'Solar', 'PMU', 'location', 'best');
% ax = gca;
% ax.FontSize = 14;
% ax.FontWeight = 'Bold';
% set(h1, 'LineWidth', 2);
% % xlim([0 1]);
% set(gcf, 'Position', [500, 300, 800, 500]);
% *********************************************************************** %

%% Plot CQI
% figure; 
% h1 = errorbar(beta, drCQI, (drCQI-drCQICI(:, 1)), (drCQICI(:, 2)-drCQI), 'b');
% grid on; hold on
% h1 = [h1; errorbar(beta, solarCQI', (solarCQI-solarCQICI(:, 1)), (solarCQICI(:, 2)-solarCQI), 'r')];
% h1 = [h1; errorbar(beta, pmuCQI', (pmuCQI-pmuCQICI(:, 1)), (pmuCQICI(:, 2)-pmuCQI), 'k')];
% xlabel('PMU \beta');
% ylabel('Channel Quality Indicator (CQI)');
% legend('DR', 'Solar', 'PMU');
% ax = gca;
% ax.FontSize = 14;
% ax.FontWeight = 'Bold';
% set(h1, 'LineWidth', 2);
% % xlim([0.1 0.9]);
% % ylim([6 9]);
% set(gcf, 'Position', [500, 300, 800, 500]);
% *********************************************************************** %

%% Delay and PDR
figure; 
yyaxis left;
h1 = errorbar(beta, pmuDelay, (pmuDelay-pmuDelayCI(:, 1)), (pmuDelayCI(:, 2)-pmuDelay));
grid on; hold on;
yyaxis right
h1 = [h1; errorbar(beta, pmuPDR, (pmuPDR-pmuPDRCI(:, 1)), (pmuPDRCI(:, 2)-pmuPDR))];
yyaxis left
ylabel('PMU average packet delay [msec]');
yyaxis right
ylabel('PMU PDR (%)');
xlabel('PMU \beta');
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'Bold';
set(h1, 'LineWidth', 1.5);
set(gcf, 'Position', [500, 300, 800, 500]);
% *********************************************************************** %

%% pmu Delay and CQI
figure; 
yyaxis left;
h1 = errorbar(beta, pmuDelay, (pmuDelay-pmuDelayCI(:, 1)), (pmuDelayCI(:, 2)-pmuDelay));
grid on; hold on;
yyaxis right
h1 = [h1; errorbar(beta, pmuCQI', (pmuCQI-pmuCQICI(:, 1)), (pmuCQICI(:, 2)-pmuCQI))];
yyaxis left
ylabel('PMU average packet delay [msec]');
yyaxis right
ylabel('PMU CQI');
xlabel('PMU \beta');
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'Bold';
set(h1, 'LineWidth', 1.5);
set(gcf, 'Position', [500, 300, 800, 500]);

%% solar Delay and CQI
figure; 
yyaxis left;
h1 = errorbar(beta, solarDelay, (solarDelay-solarDelayCI(:, 1)), (solarDelayCI(:, 2)-solarDelay));
grid on; hold on;
yyaxis right
h1 = [h1; errorbar(beta, solarCQI', (solarCQI-solarCQICI(:, 1)), (solarCQICI(:, 2)-solarCQI))];
yyaxis left
ylabel('Solar average packet delay [msec]');
yyaxis right
ylabel('solar CQI');
xlabel('PMU \beta');
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'Bold';
set(h1, 'LineWidth', 1.5);
set(gcf, 'Position', [500, 300, 800, 500]);
% *********************************************************************** %

%% PDR and CQI
figure; 
yyaxis left;
h1 = errorbar(beta, pmuPDR, (pmuPDR-pmuPDRCI(:, 1)), (pmuPDRCI(:, 2)-pmuPDR));
grid on; hold on;
yyaxis right
h1 = [h1; errorbar(beta, pmuCQI', (pmuCQI-pmuCQICI(:, 1)), (pmuCQICI(:, 2)-pmuCQI))];
yyaxis left
ylabel('PMU PDR (%)');
yyaxis right
ylabel('PMU CQI');
xlabel('PMU \beta');
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'Bold';
set(h1, 'LineWidth', 1.5);
set(gcf, 'Position', [500, 300, 800, 500]);
% *********************************************************************** %

