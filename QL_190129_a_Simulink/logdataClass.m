
% *********************************************************************** %
%% logFile class
% *********************************************************************** %

classdef logdataClass < handle
    
    properties
        
    end
    % ******************************************************************* %
    
    methods
        function obj = logdataClass
            
        end
        % *************************************************************** %
        
        %% Report function
        function report_fn(obj, varargin)
            switch (nargin-1)
                case 1
                    FileName = varargin{1};
                    
                    if(exist(FileName))
                        delete(FileName)
                    end
                case 2
                    Output      = varargin{1};
                    WhichFile   = varargin{2};
                    
                    if(~ischar(Output))
                        Output  = num2str(Output);
                        Output  = strcat(Output, ',');
                    else
                        if(Output(end) == ',')
                        else
                            if(strcmp(Output(end), 'n'))
                                % Do nothing
                            else
                                Output  = strcat(Output, ',');
                            end
                        end
                    end
                    
                    file = fopen(WhichFile, 'a');
                    fprintf(file, Output);
                    fclose(file);
                case 3
                    Output      = varargin{1};
                    WhichFile   = varargin{2};
                    format      = varargin{2};
                    
                    if(~ischar(Output))
                        Output  = num2str(Output);
                        Output  = strcat(Output, ',');
                    else
                        if(Output(end) == ',')
                        else
                            if(strcmp(Output(end), 'n'))
                                % Do nothing
                            else
                                Output  = strcat(Output, ',');
                            end
                        end
                    end
                    
                    file = fopen(WhichFile, 'a');
                    fprintf(file, format, Output);
                    fclose(file);
                otherwise
                    error('Variables for Report have errors')
            end
        end
        % *************************************************************** %

        %% Report UL Delays
        function saveLog_fn(obj, PktInfo)
            
            global logFile
            
            obj.report_fn(PktInfo.uniqueID, logFile);
            obj.report_fn(PktInfo.txType, logFile);
            obj.report_fn(PktInfo.txID, logFile);
            obj.report_fn(PktInfo.rxType, logFile);
            obj.report_fn(PktInfo.rxID, logFile);
            obj.report_fn(PktInfo.pktNumber, logFile);
            obj.report_fn(PktInfo.pktSize, logFile);
            obj.report_fn(PktInfo.crtTS, logFile);
            obj.report_fn(PktInfo.txTS, logFile);
            obj.report_fn(PktInfo.rxTS, logFile);
            obj.report_fn(PktInfo.delay, logFile);
            obj.report_fn(PktInfo.bitTp, logFile);
            obj.report_fn(PktInfo.pktErr, logFile);
            obj.report_fn('\r\n', logFile);
        end
        % *************************************************************** %
        
        %% save convergence results
        function saveconvLog_fn(obj, PktInfo)
            
            global convLogFile
            
            obj.report_fn(PktInfo.nodeID, convLogFile);
            obj.report_fn(PktInfo.nSF, convLogFile);
            obj.report_fn(PktInfo.reward, convLogFile);
            obj.report_fn(PktInfo.cqi, convLogFile);
            obj.report_fn(PktInfo.remDelay, convLogFile);
            obj.report_fn(PktInfo.prState, convLogFile);
            obj.report_fn(PktInfo.newState, convLogFile);
            obj.report_fn(PktInfo.action, convLogFile);
            obj.report_fn(PktInfo.qValue, convLogFile);
            obj.report_fn(PktInfo.epsilon, convLogFile);
            obj.report_fn('\r\n', convLogFile);
        end
        % *************************************************************** %
        
        %% save convergence results
        function storetxInfo_fn(obj, txInfo)
            
            global txlogFile
            
            obj.report_fn(txInfo.uniqueID, txlogFile);
            obj.report_fn(txInfo.nodeType, txlogFile);
            obj.report_fn(txInfo.nodeID, txlogFile);
            obj.report_fn(txInfo.crtTS, txlogFile);
            obj.report_fn(txInfo.pktNum, txlogFile);
            obj.report_fn('\r\n', txlogFile);
        end
        % *************************************************************** %
        
    end
    
end

