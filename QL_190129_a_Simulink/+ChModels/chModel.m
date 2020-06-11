
classdef chModel < handle
    
    %% DL Channel properties
    properties
        lightSpeed
        frequency
        AllFreq
        ulchcfg
        dlchcfg
        slchcfg
        delay
    end
    
    %% DL Channel methods
    methods
        %% Channel Model Constructor
        function obj = chModel()
            obj.lightSpeed = 3e8;              % 299792458
            obj.frequency = 2e9;            % 2 GHz
            obj.AllFreq = 2.5e9 + [0:49]'*180000;
            obj.delay = 25;
            
            % Channel configuration
            obj.ulchcfg.NRxAnts = 1;               % Number of receive antenna
            obj.ulchcfg.DelayProfile = 'EPA';      % Delay profile
            obj.ulchcfg.DopplerFreq = 0;           % Doppler frequency    
            obj.ulchcfg.MIMOCorrelation = 'Low';   % MIMO correlation
            obj.ulchcfg.Seed = 100;                % Channel seed    
            obj.ulchcfg.NTerms = 16;               % Oscillators used in fading model
            obj.ulchcfg.ModelType = 'GMEDS';       % Rayleigh fading model type 
            obj.ulchcfg.InitPhase = 'Random';      % Random initial phases
            obj.ulchcfg.NormalizePathGains = 'On'; % Normalize delay profile power
            obj.ulchcfg.NormalizeTxAnts = 'On';    % Normalize for transmit antennas
            
            obj.slchcfg = obj.ulchcfg;
        end
        % *************************************************************** %
        
        %% Attach Device
        function Ch_AttachDevices(obj, Device)
            Existing = find(Device == obj.Attached_Devices);
            
            if(~Existing)
                obj.Attached_Devices = [obj.Attached_Devices, Device];
            end
        end
        % *************************************************************** %
        
        %% Remove Device
        function Ch_RemoveDevices(obj, Device)
            Existing = find(Device == obj.Attached_Devices);
            
            if(Existing)
                obj.Attached_Devices(Device) = [];
            end
        end
        % *************************************************************** %
        
        %% Channel delay
        function delay = propDelay_fn(obj, distance)
            delay = (distance / obj.lightSpeed) * 1000;
        end
        % *************************************************************** %
        
        %% Apply LTE Fading channel
        function [rxWaveform, ChannelInfo] = lteCh_fn(obj, txWaveform, chType)
                        
            switch(chType)
                case 'Uplink'
                    [rxWaveform, ChannelInfo] = lteFadingChannel(obj.ulchcfg, [txWaveform; zeros(obj.delay, 1)]);
                    
                case 'Downlink'
                    [rxWaveform, ChannelInfo] = lteFadingChannel(obj.ulchcfg, [txWaveform; zeros(obj.delay, 1)]);
                    
                case 'sideLink'
                    [rxWaveform, ChannelInfo] = lteFadingChannel(obj.slchcfg, [txWaveform; zeros(obj.delay, 1)]);
            end
        end
        % *************************************************************** %
        
        %% 3GPP propagation model
        function [rxWaveform] = pathLoss_fn(obj, txWaveform, distance)
            
            % Path loss model : L = 128.1 + 37.6 log10 ( R ), R in kilometers
            % Lognormal shadowing  Log Normal Fading with 10 dB standard deviation
            
            % 3GPP pathloss (Source: Distributed Learning for Energy-Efficient
            % Resource Management in Self-Organizing Heterogeneous
            % Networks)
            distKm = distance/1000;
            PL_3GPP_dB = 128.1 + 37.6*log10(distKm);
            
            % Antenna gains
            Gt = 10;
            Gr = Gt;
            
            % Shadowing
            Mu = 0;
            Sigma_dB = 2;
            shadowingdB = normrnd(Mu, Sigma_dB);
            
            % Penetration loss and noise figure
            Penetloss_dB = 5;       % 20
            NF_dB = 5;        % Noise Figure
            
            % Total pathloss
            pathlossdB = PL_3GPP_dB + shadowingdB + Penetloss_dB + NF_dB - Gt - Gr;
            pathloss = 10^(pathlossdB/10);
            
            % Apply pathloss to the signal
            if(distKm == 0)
                rxWaveform = txWaveform;
            else
                rxWaveform = txWaveform / pathloss;
            end
        end
        % *************************************************************** %
    end
end





