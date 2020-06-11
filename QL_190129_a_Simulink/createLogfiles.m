
function createLogfiles()

%% Global variables
global logFile
global txlogFile
global convLogFile
global logData

logData = logdataClass;
% *********************************************************************** %

%% Create file to save results
logData.report_fn(logFile);

pktInfo.uniqueID = 'uniqueID';
pktInfo.txType = 'txType,';
pktInfo.txID = 'txID,';
pktInfo.rxType = 'rxType,';
pktInfo.rxID = 'rxID,';
pktInfo.pktNumber = 'pktNumber,';
pktInfo.pktSize = 'pktSize,';
pktInfo.crtTS = 'crtTS,';
pktInfo.txTS = 'txTS,';
pktInfo.rxTS = 'rxTS,';
pktInfo.delay = 'delay,';
pktInfo.bitTp = 'bitTp';    % bit throughput
pktInfo.pktErr = 'pktErr\r\n';    % bit throughput

logData.saveLog_fn(pktInfo);
% *********************************************************************** %

%% txlogFile
logData.report_fn(txlogFile);

txInfo.uniqueID = 'uniqueID,';
txInfo.nodeType = 'nodeType,';
txInfo.nodeID = 'nodeID,';
txInfo.crtTS = 'crtTS,';
txInfo.pktNum = 'pktNum\r\n';

logData.storetxInfo_fn(txInfo);
% *********************************************************************** %

%% convLogFile
logData.report_fn(convLogFile);

convInfo.nodeType = 'nodeType,';
convInfo.nodeID = 'nodeID,';
convInfo.nSF = 'nSF,';
convInfo.reward = 'reward,';
convInfo.cqi = 'cqi,';
convInfo.delay = 'delay,';
convInfo.remDelay = 'remDelay,';
convInfo.prState = 'prState,';
convInfo.newState = 'newState,';
convInfo.action = 'action,';
convInfo.qValue = 'qValue,';
convInfo.epsilon = 'epsilon\r\n';

logData.saveconvLog_fn(convInfo);
% *********************************************************************** %


end