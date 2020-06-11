
%% Test script
% This script should be the init in the simulink

clear
clc
load('HPT_1min.mat')
load('HBD_1min.mat')
load('SBD_1min.mat')
load('SPT_1min.mat')
global Seed;
Seed = 2;
rng(Seed);

%% Simulink
global activateSimulink
activateSimulink = 1;
% *********************************************************************** %

%% RunTime
if(~activateSimulink)
    global RunTime
    RunTime = 5000;      % in milli-sec
end
% *********************************************************************** %

%% Simulink globals
global rxDReventD2
global rxDReventD3
global rxDReventD4
global DRhappened
global powerReadingD2
global powerReadingD3
global powerReadingD4

rxDReventD2 = [];
rxDReventD3 = [];
rxDReventD4 = [];
DRhappened = [];
powerReadingD2 = 0;
powerReadingD3 = 0;
powerReadingD4 = 0;
% *********************************************************************** %

%% Nodes
global eNBs             % structure for eNB
global UDs              % structure for user device
global neNBs            % number of eNBs
global nD2D             % number of SGDs
global nUEs             % number of UEs

% Plot the nodes positions in the network
global networkPlot
networkPlot = 0;        

% perfectChannel
global perfectChannel
perfectChannel = 0;

% Global clock
global gClk             
gClk = 0;

% Scheduling algorithm: RoundRobin, MaxThroughput, PropFairness, QLearning,
% Random
global enbScheduler        
enbScheduler = 'Fixed'; 

global d2dScheduler
d2dScheduler = 'QLearning';
% *********************************************************************** %

%% Construct the network
Main_Init();
% *********************************************************************** %

%% Run the main
if(~activateSimulink)
    tic
    while(gClk < RunTime)
        Main();
    end
    simTime = toc;

    fprintf('\n******************************************************\n');
    fprintf('          Simulation time [min] = %d', (simTime/60));
    fprintf('\n******************************************************\n');
end
% *********************************************************************** %

%% Compute/plot the results
% global logFile 
% global convLogFile

% [solarDelay, solarCqi, solarTpt, phasorDelay, phasorCqi, ...
%     phasorTpt, drDelay, drCqi, drTpt, solarReward, phasorReward, drReward,...
%     nsfDSolar, nsfDPhasor, nsfDDr, nsfCSolar, nsfCPhasor, nsfCDr] = extractResults(logFile, convLogFile);
% 
% figure;
% scatter(nsfDSolar, solarDelay, 'b');
% title('Delay of solar panel');
% set(gcf, 'Position', [400, 200, 1000, 700]);
% 
% figure;
% scatter(nsfDPhasor, phasorDelay, 'b');
% title('Delay of phasor device');
% set(gcf, 'Position', [400, 200, 1000, 700]);
% 
% figure;
% scatter(nsfDDr, drDelay, 'b');
% title('Delay of demand-response devices');
% set(gcf, 'Position', [400, 200, 1000, 700]);

% figure;
% plot(nsfCSolar, solarCqi, 'r');
% title('CQI of solar panel');
% set(gcf, 'Position', [400, 200, 1000, 700]);
% 
% figure;
% plot(nsfCPhasor, phasorCqi, 'r');
% title('CQI of phasor device');
% set(gcf, 'Position', [400, 200, 1000, 700]);
% 
% figure;
% plot(nsfDDr, drCqi, 'r');
% title('CQI of phasor device');
% set(gcf, 'Position', [400, 200, 1000, 700]);










