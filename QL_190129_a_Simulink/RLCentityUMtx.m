
% *********************************************************************** %
% *********************************************************************** %

classdef RLCentityUMtx < handle
    
    properties
        VTus
        PDUsize
        sn_FieldLength
        
        TXbuff
        SegmentQueued
        oldest_crtTS
    end
    
    methods
       
        function obj = RLCentityUMtx(sn_FieldLength,PDUsize)
            
            if sn_FieldLength ~= 5 && sn_FieldLength ~= 10
                fprintf("WARNING: Bad sn_FieldLength input. sn_FieldLength Should be either 5 or 10. Setting sn_FieldLength to 5.\n")
                sn_FieldLength = 5;
            end
            obj.sn_FieldLength = sn_FieldLength;            
            obj.PDUsize = PDUsize;
            obj.VTus = 0;
            obj.TXbuff = {};
            obj.SegmentQueued = false;
            obj.oldest_crtTS = [];
        end
        % *************************************************************** %     
        
        function obj = rxRLCSDU(obj,sdu)
            obj.TXbuff{length(obj.TXbuff)+1} = sdu;
            if isempty(obj.oldest_crtTS)
                metadata = readMetadata(sdu,1);
                obj.oldest_crtTS = metadata.TOC;
            end
        end
        
        %%
        function [PDU,obj] = txPDU(obj,sizePDU)
            %Pumps data out of RLC TX buffer into an RLC PDU 
            %Jams as many packets as possible into RLC PDU; loosely put,
            %this involves having to at times truncate packets in buffer
            %and send what can fit into PDU while leaving the rest of it
            %behind for another PDU
            LengthLiField = 11;
            %Update TX timestamp of next in line SDU
            if ~obj.SegmentQueued
                obj.TXbuff{1} = updateTXtimestamp(obj.TXbuff{1});
                metadata = readMetadata(obj.TXbuff{1},1);
                obj.oldest_crtTS = metadata.TOC;
            end
            %Parse SDUs based on TBS: PDU size = TBS - RLC Header size
            %determine the fixed header size
            if obj.sn_FieldLength == 5
                RLCfixedHeaderSize = 8;
                frontmatter = 3;
                pduHeader = zeros(1,frontmatter);
                
                %Determine first half of FI field
                if obj.SegmentQueued
                    pduHeader(1) = 1;
                end
            else
                RLCfixedHeaderSize = 16;
                frontmatter = 6;
                pduHeader = zeros(1,frontmatter);
                
                %Determine first half of FI field
                if obj.SegmentQueued
                    pduHeader(4) = 1;
                end
            end
            
            capacity = sizePDU - RLCfixedHeaderSize;
            
            SNchar = dec2bin(obj.VTus);  
            
            %Increment sequence number
            obj.VTus = mod(obj.VTus + 1,2^obj.sn_FieldLength);
                     
            paddingLength = obj.sn_FieldLength - length(SNchar);
            SN = zeros(1,obj.sn_FieldLength);
            
            for k = paddingLength+1:obj.sn_FieldLength           
                SN(k) = str2double(SNchar(k-paddingLength));
            end
            
            pduHeader(frontmatter+1:frontmatter+length(SN)) = SN;           
            nextHeaderIndex = length(pduHeader) + 1;
            numSDUs = 0;
            nextDataIndex = 1;
            
            %Iterate through buffer until you fill a PDU or empty buffer
            while ~isempty(obj.TXbuff)
                %What is marginalHeader
                if numSDUs + 1 == 1
                    %already added fixed header
                    marginalHeader = 0;
                    lengthFirstSDU = length(obj.TXbuff{1});
                elseif numSDUs + 1 == 2
                    %weird case where you have to retroactively give a
                    %marginal header to the first SDU added in addition to
                    %marginal header for the second SDU
                    marginalHeader = 3*8;
                elseif mod(numSDUs+1,2) == 0
                    %next header will make numSDUs even
                    marginalHeader = 1.5*8;
                else
                    %next header will make numSDUs odd
                    marginalHeader = 2*8;
                end
                
                lengthSDU = length(obj.TXbuff{1});
                
                if numSDUs > 2 && mod(numSDUs+1,2) == 0
                    hypotheticalCapacity = capacity + 4;
                else
                    hypotheticalCapacity = capacity;
                end
                
                if marginalHeader + 1 < hypotheticalCapacity
                    %more data can be added to pdu
                    if numSDUs > 2 && mod(numSDUs+1,2) == 0
                        %you are about to add an SDU that makes numSDUs
                        %even, and the previous marginal header addition
                        %assumed it was the last marginal header to be
                        %added, and added 4 padding bits. Since this is not
                        %the case, these padding bits must be deleted, and
                        %nextHeaderIndex must be adjusted.
                        pduHeader(nextHeaderIndex-4:nextHeaderIndex-1) = [];
                        nextHeaderIndex = nextHeaderIndex - 4;
                        capacity = hypotheticalCapacity;
                    end
                    pduHeader(nextHeaderIndex:marginalHeader+nextHeaderIndex-1) = zeros(1,marginalHeader);
                    if lengthSDU < capacity - marginalHeader
                        %entire SDU can fit
                        capacity = capacity - lengthSDU - marginalHeader;
                        %Timestamp time of transmission
                        pduData(nextDataIndex:lengthSDU+nextDataIndex-1) = obj.TXbuff{1};
                        nextDataIndex = nextDataIndex + lengthSDU;
                        numSDUs = numSDUs + 1;
                        obj.TXbuff(1) = [];
                        obj.SegmentQueued = false;
                        %add tentative TX timestamp to next SDU in line in
                        %case it will also be put into PDU this time around
                        %also update oldest creation timestamp property
                        if ~isempty(obj.TXbuff)
                            obj.TXbuff{1} = updateTXtimestamp(obj.TXbuff{1});
                            metadata = readMetadata(obj.TXbuff{1},1);
                            obj.oldest_crtTS = metadata.TOC;
                        end
                        pduHeader(nextHeaderIndex:nextHeaderIndex+marginalHeader-1) = zeros(1,marginalHeader);
                        %populate former marginal header here?
                        
                        if numSDUs == 2 
                            %update previous header's Es AND L1 fields 
                            %Fixed header E field
                            if RLCfixedHeaderSize == 8
                                pduHeader(3) = 1;
                            else                
                                pduHeader(6) = 1;
                            end
                            %L1 E field
                            pduHeader(RLCfixedHeaderSize+1) = 1;

                            %L1 Li Field
                            %Convert lengthFirstSDU value to binary array
                            LiChar = dec2bin(lengthFirstSDU);          
                            paddingLength = LengthLiField - length(LiChar);
                            Li = zeros(1,LengthLiField);

                            for k = paddingLength+1:LengthLiField      
                                Li(k) = str2double(LiChar(k-paddingLength));
                            end

                            pduHeader(RLCfixedHeaderSize+2:RLCfixedHeaderSize+2+LengthLiField-1) = Li;
                            
                            %L2 LI Field
                            %Convert lengthSDU value to binary array
                            LiChar = dec2bin(lengthSDU);          
                            paddingLength = LengthLiField - length(LiChar);
                            Li = zeros(1,LengthLiField);

                            for k = paddingLength+1:LengthLiField      
                                Li(k) = str2double(LiChar(k-paddingLength));
                            end

                            pduHeader(RLCfixedHeaderSize+14:RLCfixedHeaderSize+14+LengthLiField-1) = Li;

                        elseif numSDUs > 2
                            %update previous header's E field
                            pduHeader(nextHeaderIndex-12) = 1;
                            
                            %update current header's LI field
                            %Convert lengthSDU value to binary array
                            LiChar = dec2bin(lengthSDU);          
                            paddingLength = LengthLiField - length(LiChar);
                            Li = zeros(1,LengthLiField);

                            for k = paddingLength+1:LengthLiField      
                                Li(k) = str2double(LiChar(k-paddingLength));
                            end

                            pduHeader(nextHeaderIndex+1:nextHeaderIndex+1+LengthLiField-1) = Li;
                        end
                        
                        nextHeaderIndex = nextHeaderIndex + marginalHeader;
                    else
                        %SDU must be segmented, and capacity is now full
                        capacity = capacity - marginalHeader;
                        pduData(nextDataIndex:capacity+nextDataIndex-1) = obj.TXbuff{1}(1:capacity);
                        obj.TXbuff{1}(1:capacity) = [];
                        obj.SegmentQueued = true;
                        numSDUs = numSDUs + 1;
                        pduHeader(nextHeaderIndex:nextHeaderIndex+marginalHeader-1) = zeros(1,marginalHeader);
                        %update FI Field in fixed header
                        if RLCfixedHeaderSize == 8
                            pduHeader(2) = 1;
                        else                
                            pduHeader(5) = 1;
                        end
                        if numSDUs == 2 
                            %update previous header's Es AND L1 fields 
                            %Fixed header E fields
                            if RLCfixedHeaderSize == 8
                                pduHeader(3) = 1;
                            else                
                                pduHeader(6) = 1;
                            end
                            %L1 E field
                            pduHeader(RLCfixedHeaderSize+1) = 1;

                            %L1 Li Field
                            %Convert lengthFirstSDU value to binary array
                            LiChar = dec2bin(lengthFirstSDU);          
                            paddingLength = LengthLiField - length(LiChar);
                            Li = zeros(1,LengthLiField);

                            for k = paddingLength+1:LengthLiField      
                                Li(k) = str2double(LiChar(k-paddingLength));
                            end

                            pduHeader(RLCfixedHeaderSize+2:RLCfixedHeaderSize+2+LengthLiField-1) = Li;
                            
                            %L2 LI Field
                            %Convert lengthSDU segment value (equals capacity) to binary array
                            LiChar = dec2bin(capacity);          
                            paddingLength = LengthLiField - length(LiChar);
                            Li = zeros(1,LengthLiField);

                            for k = paddingLength+1:LengthLiField      
                                Li(k) = str2double(LiChar(k-paddingLength));
                            end

                            pduHeader(RLCfixedHeaderSize+14:RLCfixedHeaderSize+14+LengthLiField-1) = Li;

                        elseif numSDUs > 2
                            %update previous header's E field
                            pduHeader(nextHeaderIndex-12) = 1;
                            
                            %update current header's LI field
                            %Convert lengthSDU segment value (equals capacity) to binary array
                            LiChar = dec2bin(capacity);          
                            paddingLength = LengthLiField - length(LiChar);
                            Li = zeros(1,LengthLiField);

                            for k = paddingLength+1:LengthLiField      
                                Li(k) = str2double(LiChar(k-paddingLength));
                            end

                            pduHeader(nextHeaderIndex+1:nextHeaderIndex+1+LengthLiField-1) = Li;                       
                        end
                        break
                    end
                else
                    %pdu is filled as much as possible, ready to send
                    break
                end
            end          
            PDU = [pduHeader,pduData]; 
            %after PDU is sent, check to see if tx buffer is now empty; if
            %so, delete oldest_crtTS
            if isempty(obj.TXbuff)
                obj.oldest_crtTS = [];
            end
       
            function sdu = updateTXtimestamp(sduIn)
                global gClk
                binArray = zeros(1,24);
                numStr = dec2bin(gClk);
                for i = 1:length(numStr)
                    digit = str2double(numStr(length(numStr)-i+1));
                    binArray(length(binArray)-i+1) = digit;
                end
                sduIn(1,54:77) = binArray;
                sdu = sduIn;
            end
        end
        % *************************************************************** %
        
    end
end

