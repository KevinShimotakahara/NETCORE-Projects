% *********************************************************************** %
% Dscription: Generic Packet class
% *********************************************************************** %

classdef Packet < handle
    
    properties   
        pktSize
        pktNumber
        remSize
        crtTS
        propDelay
        txTS
        pktType
    end
    
    methods
        %% Constructor
        function obj = Packet()
            obj.pktSize = [];
            obj.pktNumber = [];
            obj.remSize = [];
            obj.crtTS = [];
            obj.propDelay = [];
            obj.txTS = [];
            obj.pktType = [];
        end
        
    end
    
end


