
%% Test case:
% epsilon = 0.05 and shadowing = 2 dB
% *********************************************************************** %

%%
clear
clc

global Seed;
Seed = 2^3;
rng(Seed);
% *********************************************************************** %

%% Simulink
global activateSimulink
activateSimulink = 1;
% *********************************************************************** %

%% RunTime
global RunTime
RunTime = 7.5 * 1e3;     
% *********************************************************************** %

%% Simulink globals
if(activateSimulink)
    HPT_1min = (0:7.5/1440:7.5)';
    SPT_1min = (0:7.5/1440:7.5)';
    load('HBD_1min.mat')
    load('SBD_1min.mat')
    global DRhappened
    global powerReadingD2
    global powerReadingD3
    global powerReadingD4
    global baseComfortD2
    global baseComfortD3
    global baseComfortD4
    global ackTimer2
    global ackTimer3
    global ackTimer4
    global price
    global p
    global PriceReadingD2
    global PriceReadingD3
    global PriceReadingD4
    global demandReductionD2
    global demandReductionD3
    global demandReductionD4
    global currentDemandD2
    global currentDemandD3
    global currentDemandD4
    global reductions
    global currentDemandFunk
    
    reductions = [0,0,0];
    currentDemandD2 = [];
    currentDemandD3 = [];
    currentDemandD4 = [];
    PriceReadingD2 = [];
    PriceReadingD3 = [];
    PriceReadingD4 = [];
    ackTimer2 = [];
    ackTimer3 = [];
    ackTimer4 = [];
    baseComfortD2 = 320;
    baseComfortD3 = 320;
    baseComfortD4 = 320;
    DRhappened = [];
    powerReadingD2 = 0;
    powerReadingD3 = 0;
    powerReadingD4 = 0;
    demandReductionD2 = 0;
    demandReductionD3 = 0;
    demandReductionD4 = 0;
    
    time = (0:2:7500)';
%     priceMax = 500;
%     priceMin = 200;
%     tmed = 3750;
%     p = @(t) ((priceMin-priceMax)/(tmed^2))*(t-tmed)^2 + priceMax;
    for k = 1:length(time)
        if time(k) <= 2000
            priceVec(k) = 200;
        elseif time(k) <= 2500
            priceVec(k) = 300;
        elseif time(k) <= 3000
            priceVec(k) = 500;
        else
            priceVec(k) = 200;
        end
    end
    priceVec = priceVec';
    price = [time/1000 priceVec];
    p = @(t) interp1(time,priceVec,t);
    currentDemandFunk = @(t) interp1(HPT_1min,HPD_1min,t);
end
% *********************************************************************** %

%% Nodes
global eNBs             % structure for eNB
global UDs              % structure for user device
global neNBs            % number of eNBs
global nD2D             % number of SGDs
global nUEs             % number of UEs

% perfectChannel
global perfectChannel
perfectChannel = 0;

% Global clock
global gClk             
gClk = 0;

% Run number
global nRun
nRun = 1;

% Scheduling algorithm: RoundRobin, MaxThroughput, PropFairness, QLearning,
% Random
global enbScheduler        
enbScheduler = 'Fixed'; 

global d2dScheduler
d2dScheduler = 'RandomNew';
% *********************************************************************** %

%% Construct the network
Main_Init();
% *********************************************************************** %

%% Draw the positions
global networkPlot
networkPlot = 0;   
PlotNetPlan(networkPlot);
% *********************************************************************** %

%% Run the main
if(~activateSimulink)
    tic
    while(nRun <= 15)
        Main();
    end
    
    simTime = (toc/3600);
    fprintf('\n******************************************************\n');
    fprintf('          Simulation time [hours] = %d', simTime);
    fprintf('\n******************************************************\n');
    
    save('simTime', 'simTime');
end
% *********************************************************************** %


