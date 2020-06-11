
function Init_simParameters()

global simParameters

%% General parameters
simParameters.bandwidth = 5;           % in MHz
simParameters.RBBW = 180000;            % Resource block bandwidth (Hz)
simParameters.nRBs = 50;

simParameters.enbMinSeparation = 600;   % in meters
simParameters.cellRadius = 500;         % simParameters.enbSeparation / 2;

simParameters.pathloss = '3GPP';        % Pathloss model

% Traffic pattern
ipHeader = 77;
simParameters.ueTraffic.avgPktSize = 2000-ipHeader; 
simParameters.ueTraffic.lambda = 20;            %1/40e-3: 1 packets/ 40 milli-sec
simParameters.drTraffic.avgPktSize = 2000-ipHeader;
simParameters.drTraffic.lambda = 20;       
simParameters.solarTraffic.avgPktSize = 1120-ipHeader;
simParameters.solarTraffic.lambda = 20;        
simParameters.pmuTraffic.avgPktSize = 592-ipHeader;
simParameters.pmuTraffic.lambda = 16;       

% Delay requirement
simParameters.dummyTraffic.delayMin = 1200;
simParameters.dummyTraffic.delayMax = 1000;
simParameters.drTraffic.delayMin = 500;
simParameters.drTraffic.delayMax = 600;
simParameters.solarTraffic.delayMin = 300;
simParameters.solarTraffic.delayMax = 500;
simParameters.pmuTraffic.delayMin = 30;
simParameters.pmuTraffic.delayMax = 50;

% Beta configuration
simParameters.dummyTraffic.beta = 0;
simParameters.drTraffic.beta = 0;         
simParameters.solarTraffic.beta = 0;
simParameters.pmuTraffic.beta = 0;
% *********************************************************************** %

%% ulConfig parameters (UE)
ue.NCellID = 10;     % Cell identity
ue.RC = 'A3-7';      % FRC number

% UE ulConfig: FRC configuration structure A3-7 (Just an init)
frc = lteRMCUL(ue);
ulConfig = frc;
ulConfig.TotSubframes = 1;

% ue HARQ init
ulConfig.HARQ.harqProcessSequence = 1;
ulConfig.HARQ.harqProcesses = hPUSCHNewHARQProcess(ulConfig); 
ulConfig.HARQ.harqID = 1;
ulConfig.HARQ.harqProcIdx = 0;
ulConfig.sdus = {};

simParameters.ue = ue;
simParameters.ulConfig = ulConfig;
% *********************************************************************** %

%% dlConfig parameters (eNB)

% *********************************************************************** %

%% slConfig parameters (D2D)
simParameters.slConfig = ulConfig;
simParameters.slConfig.sciMessageRx = [];
simParameters.slConfig.sciDecoded = 0;
simParameters.slConfig.cqiBuffer = [];
simParameters.slConfig.recentDelay = 100;
% *********************************************************************** %

end




