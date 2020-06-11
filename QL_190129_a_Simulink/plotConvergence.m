

%% Variables
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

beta = 0;
% *********************************************************************** %

%% Convergence for one run
[nSFdr1, nSFdr2, nSFdr3, nSFdr4, nSFsolar, nSFpmu,...
    dr1Data, dr2Data, dr3Data, dr4Data, solarData, pmuData] = extractBetaConv(6, beta);

% DR 1
figure; 
h1 = scatter(nSFdr1, dr1Data, 'b');
grid on;
ylabel('DR 1 packet delay [msec]');
xlabel('Subframe number');
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'Bold';
set(h1, 'LineWidth', 1.5);
set(gcf, 'Position', [500, 300, 800, 500]);

% DR 2
figure; 
h1 = scatter(nSFdr2, dr2Data, 'b');
grid on;
ylabel('DR 2 packet delay [msec]');
xlabel('Subframe number');
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'Bold';
set(h1, 'LineWidth', 1.5);
set(gcf, 'Position', [500, 300, 800, 500]);

% DR 3
figure; 
h1 = scatter(nSFdr3, dr3Data, 'b');
grid on;
ylabel('DR 3 packet delay [msec]');
xlabel('Subframe number');
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'Bold';
set(h1, 'LineWidth', 1.5);
set(gcf, 'Position', [500, 300, 800, 500]);

% DR 4
figure; 
h1 = scatter(nSFdr4, dr4Data, 'b');
grid on;
ylabel('DR 4 packet delay [msec]');
xlabel('Subframe number');
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'Bold';
set(h1, 'LineWidth', 1.5);
set(gcf, 'Position', [500, 300, 800, 500]);

% Solar
figure; 
h1 = scatter(nSFsolar, solarData, 'b');
grid on;
ylabel('Solar packet delay [msec]');
xlabel('Subframe number');
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'Bold';
set(h1, 'LineWidth', 1.5);
set(gcf, 'Position', [500, 300, 800, 500]);

% PMU
figure; 
h1 = scatter(nSFpmu, pmuData, 'b');
grid on;
ylabel('PMU packet delay [msec]');
xlabel('Subframe number');
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'Bold';
set(h1, 'LineWidth', 1.5);
set(gcf, 'Position', [500, 300, 800, 500]);
% *********************************************************************** %









