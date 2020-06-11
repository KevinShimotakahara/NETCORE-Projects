
% Initlization function for the global variable.

function Main_Init()

%% Global variables
global nRun

global neNBs        % Number of eNBs in the network.
global nD2D
global nDummy
global nUEs
global nUDs

global UDs          % Node structure for UDs
global eNBs         % Node structure for eNB

global simParameters
global enbScheduler
global d2dScheduler
global dlChannel
global ulChannel
% *********************************************************************** %

%% Initialize the simulation
simParameters  = struct;
Init_simParameters();
% *********************************************************************** %

%% Construct network elements
neNBs = 1;          % Number of eNBs
nD2D = 7;           % Number of D2D users per eNB
nDummy = 5;
nUEs = 0;           % Number of UEs users per eNB
nUDs = nD2D + (2*nDummy); % total number of users per eNB

UDs = NetElements.UD;       % struct for users
eNBs = NetElements.eNodeB;  % struct for eNBs
% *********************************************************************** %

%% Log file: ".dat" file to save the results  
global logFile      
logFile = strcat('./Results/', num2str(neNBs), 'eNBs_', num2str(nD2D),...
'D2Ds_', num2str(nDummy), 'Dummy_', num2str(nUEs), 'UEs_', d2dScheduler, '_Run', num2str(nRun),...
'_Beta', num2str(simParameters.pmuTraffic.beta), '.dat');

global txlogFile      
txlogFile = strcat('./Results/tx_', num2str(neNBs), 'eNBs_', num2str(nD2D),...
'D2Ds_', num2str(nDummy), 'Dummy_', num2str(nUEs), 'UEs_', d2dScheduler, '_Run', num2str(nRun),...
'_Beta', num2str(simParameters.pmuTraffic.beta), '.dat');

global convLogFile  % A file for saving the convergence results
convLogFile = strcat('./Results/convLogFile_', num2str(neNBs), 'eNBs_', num2str(nD2D),...
'D2Ds_', num2str(nDummy), 'Dummy_', num2str(nUEs), 'UEs_', d2dScheduler, '_Run', num2str(nRun),...
'_Beta', num2str(simParameters.pmuTraffic.beta), '.dat');
% *********************************************************************** %

%% Construct and initialize the eNBs
for k = 1:neNBs
    eNBs(k).nodeID = k;
    eNBs(k).nodeType = 'eNB';
    eNBs(k).cellRadius = 50; % 1 km: Power control in two-tier femtocell networks

    eNBs(k).setScheduler_fn(enbScheduler);         % Schedulers: RoundRobin, QLearning, PropFairness, MaxThroughput
    eNBs(k).enbScheduler.agentNode = eNBs(k);   % the node running the scheduler
end

simParameters.GridSize = 4*max([eNBs.cellRadius])/1000;
% *********************************************************************** %

%% Create users (D2D, and UEs) 
for u = 1:nUDs
    UDs(u) = NetElements.UD;
end

% Allocate txUsers and rxUsers
UDs(1).rxUsers = [UDs(2); UDs(3); UDs(4)];      % users receiving from 1 (2, 3, and 4 transmit to 1)
UDs(2).txUsers = UDs(1);        % user 1 transmits to 2
UDs(3).txUsers = UDs(1);        % user 1 transmits to 3
UDs(4).txUsers = UDs(1);        % user 1 transmits to 4

UDs(1).txUsers = [UDs(2); UDs(3); UDs(4)];      % users transmitting to 1 (2, 3, and 4 transmit to 1)
UDs(2).rxUsers = UDs(1);        % user 2 receives from 1
UDs(3).rxUsers = UDs(1);        % user 3 receives from 1
UDs(4).rxUsers = UDs(1);        % user 4 receives from 1

UDs(5).txUsers = UDs(7);        
UDs(7).rxUsers = UDs(5);

UDs(6).txUsers = UDs(5);        % user 5 transmits to user 6
UDs(5).rxUsers = UDs(6);       

% Configure users for transmission
% UDs(1).txDisable = 1;
% UDs(2).txDisable = 1;
% UDs(3).txDisable = 1;
% UDs(4).txDisable = 1;
% UDs(5).txDisable = 0;
% UDs(6).txDisable = 0;
% UDs(7).txDisable = 0;

% Pool configuration (fixed)
nRBGs = 5;
nStates = 6;

i = 1;
udIdx = 1;

% (ID: 1) transformer: tx/rx to Demand-Response (DR) user
UDs(udIdx).nodeID = udIdx;
UDs(udIdx).absID = udIdx;
UDs(udIdx).nodeType = 'D2D';      
UDs(udIdx).D2DType = 'DR';      
UDs(udIdx).commDir = 'txrx';  
UDs(udIdx).cellID = i;
UDs(udIdx).cellType = 'eNB';
UDs(udIdx).cellAttached = eNBs(i);      % base station that covers this user.
UDs(udIdx).txPower = 20;       % tx Power in dB
UDs(udIdx).ulConfig = [];      % uplink configuration (no uplink for D2D)
UDs(udIdx).dlConfig = [];      % downlink configuration (no downlink for D2D)
UDs(udIdx).sltxConfig = simParameters.slConfig;         % sidelink configuration 
UDs(udIdx).RLCtxEntity = RLCentityUMtx(10, 772); % RLCentityUMtx(sn_FieldLength, sizePDU);
for u = 1:3     % receive entity for each Demand-Response user
    UDs(udIdx).RLCrxEntity{u} = RLCentityUMrx(10, 772, 1); % RLCentityUMtx(sn_FieldLength, sizePDU);
    UDs(udIdx).slrxConfig{u} = simParameters.slConfig;         % sidelink configuration 
end
UDs(udIdx).pscchPeriod = pscch;
UDs(udIdx).allTraffic_fn(simParameters.drTraffic.lambda); % Generate traffic for the whole simulation steps.    
UDs(udIdx).setScheduler_fn(d2dScheduler);
UDs(udIdx).d2dScheduler.initScheduler_fn(nRBGs, nStates, simParameters.drTraffic.delayMin,...
    simParameters.drTraffic.delayMax, 9, 9, UDs(udIdx), simParameters.drTraffic.beta);   
eNBs(i).udsAttach_fn(UDs(udIdx));

udIdx = udIdx + 1;

% (ID: 2, 3, 4) Demand-Response users: tx/rx to transformer 1
for j = 1:3
    UDs(udIdx).nodeID = udIdx;
    UDs(udIdx).absID = j;
    UDs(udIdx).nodeType = 'D2D';  
    UDs(udIdx).D2DType = 'DR';      
    UDs(udIdx).commDir = 'txrx';   
    UDs(udIdx).cellID = i;
    UDs(udIdx).cellType = 'eNB';
    UDs(udIdx).cellAttached = eNBs(i);      % base station that covers this user.
    UDs(udIdx).txPower = 20;       % tx Power in dB
    UDs(udIdx).ulConfig = [];      % uplink configuration (no uplink for D2D)
    UDs(udIdx).dlConfig = [];      % downlink configuration (no downlink for D2D)
    UDs(udIdx).sltxConfig = simParameters.slConfig;         % sidelink configuration 
    UDs(udIdx).RLCtxEntity = RLCentityUMtx(10, 772); % RLCentityUMtx(sn_FieldLength, sizePDU);
    UDs(udIdx).RLCrxEntity{1} = RLCentityUMrx(10, 772, 1); % RLCentityUMtx(sn_FieldLength, sizePDU);
    UDs(udIdx).slrxConfig{1} = simParameters.slConfig;         % sidelink configuration 
    UDs(udIdx).pscchPeriod = pscch;
    UDs(udIdx).allTraffic_fn(simParameters.drTraffic.lambda); % Generate traffic for the whole simulation steps.    
    UDs(udIdx).setScheduler_fn(d2dScheduler);
    UDs(udIdx).d2dScheduler.initScheduler_fn(nRBGs, nStates, simParameters.drTraffic.delayMin,...
        simParameters.drTraffic.delayMax, 9, 9, UDs(udIdx), simParameters.drTraffic.beta);   
    eNBs(i).udsAttach_fn(UDs(udIdx));
    
    udIdx = udIdx + 1;
end

% (ID: 5) PMU: tx to node 6 (PWM), rx from 7 (solar)
UDs(udIdx).nodeID = udIdx;
UDs(udIdx).absID = 1;
UDs(udIdx).nodeType = 'D2D';  
UDs(udIdx).D2DType = 'PMU';      
UDs(udIdx).commDir = 'txrx';  
UDs(udIdx).cellID = i;
UDs(udIdx).cellType = 'eNB';
UDs(udIdx).cellAttached = eNBs(i);      % base station that covers this user.
UDs(udIdx).txPower = 20;       % tx Power in dB
UDs(udIdx).ulConfig = [];      % uplink configuration (no uplink for D2D)
UDs(udIdx).dlConfig = [];      % downlink configuration (no downlink for D2D)
UDs(udIdx).sltxConfig = simParameters.slConfig;         % sidelink configuration 
UDs(udIdx).RLCtxEntity = RLCentityUMtx(10, 772); % RLCentityUMtx(sn_FieldLength, sizePDU);
UDs(udIdx).RLCrxEntity{1} = RLCentityUMrx(10, 772, 1); % RLCentityUMtx(sn_FieldLength, sizePDU);
UDs(udIdx).slrxConfig{1} = simParameters.slConfig;         % sidelink configuration 
UDs(udIdx).pscchPeriod = pscch;
UDs(udIdx).allTraffic_fn(simParameters.pmuTraffic.lambda); % Generate traffic for the whole simulation steps.    
UDs(udIdx).setScheduler_fn(d2dScheduler);
UDs(udIdx).d2dScheduler.initScheduler_fn(nRBGs, nStates, simParameters.pmuTraffic.delayMin,...
    simParameters.pmuTraffic.delayMax, 9, 9, UDs(udIdx), simParameters.pmuTraffic.beta);   
eNBs(i).udsAttach_fn(UDs(udIdx));

udIdx = udIdx + 1;

% (ID: 6) PWM controller: rx from PMU
UDs(udIdx).nodeID = udIdx;
UDs(udIdx).absID = 1;
UDs(udIdx).nodeType = 'D2D'; 
UDs(udIdx).commDir = 'rx';   
UDs(udIdx).cellID = i;
UDs(udIdx).cellType = 'eNB';
UDs(udIdx).cellAttached = eNBs(i);      % base station that covers this user.
UDs(udIdx).txPower = 20;       % tx Power in dB
UDs(udIdx).ulConfig = [];      % uplink configuration (no uplink for D2D)
UDs(udIdx).dlConfig = [];      % downlink configuration (no downlink for D2D)
UDs(udIdx).sltxConfig = simParameters.slConfig;         % sidelink configuration 
UDs(udIdx).RLCtxEntity = RLCentityUMtx(10, 772); % RLCentityUMtx(sn_FieldLength, sizePDU);
UDs(udIdx).RLCrxEntity{1} = RLCentityUMrx(10, 772, 1); % RLCentityUMtx(sn_FieldLength, sizePDU);
UDs(udIdx).slrxConfig{1} = simParameters.slConfig;         % sidelink configuration 
UDs(udIdx).pscchPeriod = pscch;
UDs(udIdx).allTraffic_fn(simParameters.pmuTraffic.lambda); % Generate traffic for the whole simulation steps.    
UDs(udIdx).setScheduler_fn(d2dScheduler);
UDs(udIdx).d2dScheduler.initScheduler_fn(nRBGs, nStates, simParameters.solarTraffic.delayMin,...
    simParameters.solarTraffic.delayMax, 9, 9, UDs(udIdx), simParameters.dummyTraffic.beta); 
eNBs(i).udsAttach_fn(UDs(udIdx));

udIdx = udIdx + 1;

% (ID: 7) Solar Panel: tx to node 5 (PMU)
UDs(udIdx).nodeID = udIdx;
UDs(udIdx).absID = 1;
UDs(udIdx).nodeType = 'D2D'; 
UDs(udIdx).D2DType = 'Solar';      
UDs(udIdx).commDir = 'tx';    
UDs(udIdx).cellID = i;
UDs(udIdx).cellType = 'eNB';
UDs(udIdx).cellAttached = eNBs(i);      % base station that covers this user.
UDs(udIdx).txPower = 20;       % tx Power in dB
UDs(udIdx).ulConfig = [];      % uplink configuration (no uplink for D2D)
UDs(udIdx).dlConfig = [];      % downlink configuration (no downlink for D2D)
UDs(udIdx).sltxConfig = simParameters.slConfig;         % sidelink configuration 
UDs(udIdx).RLCtxEntity = RLCentityUMtx(10, 772); % RLCentityUMtx(sn_FieldLength, sizePDU);
UDs(udIdx).RLCrxEntity{1} = RLCentityUMrx(10, 772, 1); % RLCentityUMtx(sn_FieldLength, sizePDU);
UDs(udIdx).slrxConfig{1} = simParameters.slConfig;         % sidelink configuration 
UDs(udIdx).pscchPeriod = pscch;
UDs(udIdx).allTraffic_fn(simParameters.solarTraffic.lambda); % Generate traffic for the whole simulation steps.    
UDs(udIdx).setScheduler_fn(d2dScheduler);
UDs(udIdx).d2dScheduler.initScheduler_fn(nRBGs, nStates, simParameters.solarTraffic.delayMin,...
    simParameters.solarTraffic.delayMax, 9, 9, UDs(udIdx), simParameters.solarTraffic.beta);   % init the scheduler
eNBs(i).udsAttach_fn(UDs(udIdx));

udIdx = udIdx + 1;

%% Dummy devices
for j = 1:nDummy
    UDs(udIdx).nodeID = udIdx;
    UDs(udIdx).absID = 1;
    UDs(udIdx).nodeType = 'DD2D';
    UDs(udIdx).nodeDummy = nDummy;
    UDs(udIdx).commDir = 'tx';     
    UDs(udIdx).cellID = i;
    UDs(udIdx).cellType = 'eNB';
    UDs(udIdx).cellAttached = eNBs(i);      % base station that covers this user.
    UDs(udIdx).txPower = 20;       % tx Power in dB
    UDs(udIdx).ulConfig = [];      % uplink configuration (no uplink for D2D)
    UDs(udIdx).dlConfig = [];      % downlink configuration (no downlink for D2D)
    UDs(udIdx).sltxConfig = simParameters.slConfig;         % sidelink configuration 
    UDs(udIdx).RLCtxEntity = RLCentityUMtx(10, 772); % RLCentityUMtx(sn_FieldLength, sizePDU);
    UDs(udIdx).RLCrxEntity{1} = RLCentityUMrx(10, 772, 1); % RLCentityUMtx(sn_FieldLength, sizePDU);
    UDs(udIdx).slrxConfig{1} = simParameters.slConfig;         % sidelink configuration 
    UDs(udIdx).pscchPeriod = pscch;
    UDs(udIdx).allTraffic_fn(simParameters.ueTraffic.lambda); % Generate traffic for the whole simulation steps.    
    UDs(udIdx).setScheduler_fn(d2dScheduler);
    UDs(udIdx).d2dScheduler.initScheduler_fn(nRBGs, nStates, simParameters.dummyTraffic.delayMin,...
    simParameters.dummyTraffic.delayMax, 9, 9, UDs(udIdx), simParameters.dummyTraffic.beta);   % init the scheduler
    eNBs(i).udsAttach_fn(UDs(udIdx));
    
    UDs(udIdx).rxUsers = UDs(udIdx+1);

    udIdx = udIdx + 1;
    
    UDs(udIdx).nodeID = udIdx;
    UDs(udIdx).absID = 1;
    UDs(udIdx).nodeType = 'DD2D'; 
    UDs(udIdx).nodeDummy = nDummy;
    UDs(udIdx).commDir = 'rx';     
    UDs(udIdx).cellID = i;
    UDs(udIdx).cellType = 'eNB';
    UDs(udIdx).cellAttached = eNBs(i);      % base station that covers this user.
    UDs(udIdx).txPower = 20;       % tx Power in dB
    UDs(udIdx).ulConfig = [];      % uplink configuration (no uplink for D2D)
    UDs(udIdx).dlConfig = [];      % downlink configuration (no downlink for D2D)
    UDs(udIdx).sltxConfig = simParameters.slConfig;         % sidelink configuration 
    UDs(udIdx).RLCtxEntity = RLCentityUMtx(10, 772); % RLCentityUMtx(sn_FieldLength, sizePDU);
    UDs(udIdx).RLCrxEntity{1} = RLCentityUMrx(10, 772, 1); % RLCentityUMtx(sn_FieldLength, sizePDU);
    UDs(udIdx).slrxConfig{1} = simParameters.slConfig;         % sidelink configuration 
    UDs(udIdx).pscchPeriod = pscch;
    UDs(udIdx).allTraffic_fn(simParameters.ueTraffic.lambda); % Generate traffic for the whole simulation steps.    
    UDs(udIdx).setScheduler_fn(d2dScheduler);
    UDs(udIdx).d2dScheduler.initScheduler_fn(nRBGs, nStates, simParameters.dummyTraffic.delayMin,...
    simParameters.dummyTraffic.delayMax, 9, 9, UDs(udIdx), simParameters.dummyTraffic.beta);   % init the scheduler
    eNBs(i).udsAttach_fn(UDs(udIdx));
    
    UDs(udIdx).txUsers = UDs(udIdx-1);

    udIdx = udIdx + 1;
end
% *********************************************************************** %

%% Allocate positions to network elements

% Allocate eNBs positions
for e = 1:length(eNBs)
    eNBs(e).enbPosalloc_fn();
end

% Allocate UDs positions
for e = 1 : length(eNBs)
    eNBs(e).udsPosalloc_fn();
end
% ********************************************************************** %

%% Update interferers list of the receiving nodes
UDs(1).udsInterfering = [UDs(1), UDs(2), UDs(3), UDs(4), UDs(5), UDs(7)]; % interferers list of receiving node 1
UDs(2).udsInterfering = [UDs(1), UDs(2), UDs(3), UDs(4), UDs(5), UDs(7)];
UDs(3).udsInterfering = [UDs(1), UDs(2), UDs(3), UDs(4), UDs(5), UDs(7)];
UDs(4).udsInterfering = [UDs(1), UDs(2), UDs(3), UDs(4), UDs(5), UDs(7)];
UDs(5).udsInterfering = [UDs(1), UDs(2), UDs(3), UDs(4), UDs(5)];
UDs(6).udsInterfering = [UDs(1), UDs(2), UDs(3), UDs(4), UDs(7)];

for u = 1:length(UDs)
    if( strcmp(UDs(u).nodeType, 'DD2D') && (strcmp(UDs(u).commDir, 'tx') || strcmp(UDs(u).commDir, 'txrx')) )
        UDs(1).udsInterfering = [UDs(1).udsInterfering, UDs(u)];
        UDs(2).udsInterfering = [UDs(2).udsInterfering, UDs(u)];
        UDs(3).udsInterfering = [UDs(3).udsInterfering, UDs(u)];
        UDs(4).udsInterfering = [UDs(4).udsInterfering, UDs(u)];
        UDs(5).udsInterfering = [UDs(5).udsInterfering, UDs(u)];
        UDs(6).udsInterfering = [UDs(6).udsInterfering, UDs(u)];
    end
    
    if( strcmp(UDs(u).nodeType, 'DD2D') && (strcmp(UDs(u).commDir, 'rx') || strcmp(UDs(u).commDir, 'txrx')) ) 
        otherDummies = find( (strcmp({UDs.nodeType}, 'DD2D') & strcmp({UDs.commDir}, 'tx')) );
        otherDummies(otherDummies == u) = [];
        otherDummies(otherDummies == u-1) = []; %MAY NEED TO CHANGE THIS IF DD2D NETWORK CHANGES
        UDs(u).udsInterfering = [UDs(1), UDs(2), UDs(3), UDs(4), UDs(5), UDs(7), UDs(otherDummies)];
    end
end
% *********************************************************************** %

%% Update distance to D2D rxer
for u = 1:length(UDs)
    udTemp = UDs(u);
    if(strcmp(udTemp.nodeType, 'D2D') || strcmp(udTemp.nodeType, 'DD2D'))
        udTemp.disttoRxer_fn(udTemp.rxUsers);  
    end
end
% *********************************************************************** %

%% Channel model
dlChannel = ChModels.chModel;
ulChannel = ChModels.chModel;
% *********************************************************************** %

%% Create log files
createLogfiles();
% *********************************************************************** %

end











