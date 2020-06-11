
% *********************************************************************** %
% *********************************************************************** %

classdef RLCentityUMrx < handle

    properties
        VTus %This state variable holds the value of the SN to be assigned for the next newly generated UMD PDU. It is initially set to
             %0, and is updated whenever the UM RLC entity delivers an UMD PDU with SN = VT(US). 

        VRur %This state variable holds the value of the SN of the earliest UMD PDU that is still considered for reordering. It is
             %initially set to 0.

        VRuh %This state variable holds the value of the SN following the SN of the UMD PDU with the highest SN among received
             %UMD PDUs, and it serves as the higher edge of the reordering window. It is initially set to 0. 

        VRux %This state variable holds the value of the SN following the SN of the UMD PDU 

        UM_Window_Size  %This constant is used by the receiving UM RLC entity to define SNs of those UMD PDUs that can be received without
                        %causing an advancement of the receiving window. UM_Window_Size = 16 when a 5 bit SN is configured, 
                        %UM_Window_Size = 512 when a 10 bit SN is configured and UM_Window_Size = 0 when the receiving UM RLC
                        %entity is configured for MCCH, MTCH, SC-MCCH, SC-MTCH or STCH.

        T_reordering %This timer is used by the receiving side of an AM RLC entity and receiving UM RLC entity in order to detect loss of
                     %RLC PDUs at lower layer (see sub clauses 5.1.2.2 and 5.1.3.2). If t-Reordering is running, t-Reordering shall not be
                     %started additionally, i.e. only one t-Reordering per RLC entity is running at a given time. 
        PDUsize
        sn_FieldLength

        SDUs %Reassembled RLC SDUs are temporarily stored in SDUs property so they can be accessed by higher layer protocols via the RLC rx entity object
        reassemblyOccured %allows external protocols to know when new SDUs were unpacked by RLC layer
        RXbuff
        SegmentQueued
        someresetvalue
        SDUfragment
        SDUfragmentSN
        SDUsDropped
    end

    methods
       
        function obj = RLCentityUMrx(sn_FieldLength,PDUsize,someresetvalue)
            if sn_FieldLength == 5
                obj.UM_Window_Size = 16;
            elseif sn_FieldLength == 10
                obj.UM_Window_Size = 512;
            else
                fprintf("WARNING: Bad sn_FieldLength input. sn_FieldLength Should be either 5 or 10. Setting sn_FieldLength to 5.\n")
                sn_FieldLength = 5;
                obj.UM_Window_Size = 16;
            end
            
            obj.VRuh = 0;
            obj.PDUsize = PDUsize;
            obj.sn_FieldLength = sn_FieldLength;
            obj.VRur = 0;
            obj.VTus = 0;
            obj.VRux = 'NULL'; %This state variable holds the value of the SN following the SN of the UMD PDU 
            obj.SegmentQueued = false;
            obj.T_reordering.running = false;
            obj.RXbuff = {};
            obj.someresetvalue = someresetvalue;
            obj.T_reordering.timer = someresetvalue;
            obj.SDUs = {};
            obj.reassemblyOccured = false;
            obj.SDUfragment = [];
            obj.SDUfragmentSN = -1;
            obj.SDUsDropped = 0;
        end
        % *************************************************************** %     
        
        %%
        function obj = rxPDU(obj,PDU)
            obj.SDUsDropped = 0;
            PDU = double(PDU);
            %modulo thing
            moduloDivisor = 2^obj.sn_FieldLength;
            %Need to see if you will even accept incoming PDU
            x = readRLCPDUSN(obj,PDU);
            rxBefore = false;
            %Ignore PDU with same SN as one in RXbuff (already have this
            %PDU)
            for k = 1:length(obj.RXbuff)
                if modulothing(x,obj,moduloDivisor) == modulothing(readRLCPDUSN(obj,obj.RXbuff{k}),obj,moduloDivisor)
                    rxBefore = true;
                    break
                end
            end
            
            if ((modulothing(obj.VRur,obj,moduloDivisor) < modulothing(x,obj,moduloDivisor)) && (modulothing(x,obj,moduloDivisor)...
                    < modulothing(obj.VRuh,obj,moduloDivisor)) && rxBefore) ||  ...
                    ((0 <= modulothing(x,obj,moduloDivisor)) && (modulothing(x,obj,moduloDivisor) < modulothing(obj.VRur,obj,moduloDivisor)))
                return; %packet is discarded
            else
                %add to buffer
                obj.RXbuff{length(obj.RXbuff)+1} = PDU;
            end
            
            if modulothing(x,obj,moduloDivisor) >= modulothing(obj.VRuh,obj,moduloDivisor)
                %x falls outside reordering window
                obj.VRuh = mod(x + 1,moduloDivisor);
                %reassemble RLC SDUs from any UMD PDUs with SN that falls outside of the reordering window, remove
                %RLC headers when doing so and deliver the reassembled RLC SDUs to upper layer in ascending order of the
                %RLC SN if not delivered before; 
                obj = reassemble(obj,2);            
                if modulothing(obj.VRur,obj,moduloDivisor) >= modulothing(obj.VRuh,obj,moduloDivisor)
                    %VR(UR) falls outside of the reordering window:
                    obj.VRur = mod(obj.VRuh - obj.UM_Window_Size,moduloDivisor);
                end
            end      
            
            %See if the reception buffer contains an UMD PDU with SN = VR(UR): 
            VRurRXed = false;
            RXbuffSNs = zeros(length(obj.RXbuff),2);
            for k = 1:length(obj.RXbuff)
                RXbuffSNs(k,2) = readRLCPDUSN(obj,obj.RXbuff{k});
                RXbuffSNs(k,1) = modulothing(RXbuffSNs(k,2),obj,moduloDivisor);
                if modulothing(obj.VRur,obj,moduloDivisor) == RXbuffSNs(k,1)
                    VRurRXed = true;
                end
            end
            
            latch = false;
            if VRurRXed
                %update VR(UR) to the SN of the first UMD PDU with SN > current VR(UR) that has not been received
                if max(RXbuffSNs(:,1)) == modulothing(obj.VRur,obj,moduloDivisor)
                    obj.VRur = mod(obj.VRur + 1,moduloDivisor);
                else
                    RXbuffSNs = sortrows(RXbuffSNs);
                    if length(RXbuffSNs(:,1)) == 1
                        lastk = 1;
                    else
                        lastk = length(RXbuffSNs(:,1))-1;
                    end
                    for k = 1:lastk
                        if RXbuffSNs(k,1) > modulothing(obj.VRur,obj,moduloDivisor)
                            if latch == false
                                latch = true;
                                if RXbuffSNs(k,1)-1 ~= modulothing(obj.VRur,obj,moduloDivisor)
                                    obj.VRur = mod(obj.VRur + 1,moduloDivisor);
                                    break
                                elseif RXbuffSNs(k,1)+1 ~= RXbuffSNs(k+1,1)
                                    obj.VRur = mod(RXbuffSNs(k,2)+1,moduloDivisor);
                                    break
                                elseif k == length(RXbuffSNs) - 1
                                    obj.VRur = mod(RXbuffSNs(k+1,2)+1,moduloDivisor);
                                end
                            end
                            if RXbuffSNs(k,1)+1 ~= RXbuffSNs(k+1,1)
                                obj.VRur = mod(RXbuffSNs(k,2)+1,moduloDivisor);
                                break
                            elseif k == length(RXbuffSNs) - 1
                                obj.VRur = mod(RXbuffSNs(k+1,2)+1,moduloDivisor);
                            end
                        end
                    end
                end
                obj = reassemble(obj,1);
            end

            if obj.T_reordering.running
                % - if VR(UX) <= VR(UR); or
                % - if VR(UX) falls outside of the reordering window and VR(UX) is not equal to VR(UH):: 
                if modulothing(obj.VRux,obj,moduloDivisor) <= modulothing(obj.VRur,obj,moduloDivisor) ||...
                        (modulothing(obj.VRux,obj,moduloDivisor) >= modulothing(obj.VRuh,obj,moduloDivisor) ...
                        && modulothing(obj.VRux,obj,moduloDivisor) ~= modulothing(obj.VRuh,obj,moduloDivisor))
                    %stop and reset t-Reordering; 
                    obj.T_reordering.running = false;
                    obj.T_reordering.timer = obj.someresetvalue;
                end
            end
            if obj.T_reordering.running == false
                if modulothing(obj.VRuh,obj,moduloDivisor) > modulothing(obj.VRur,obj,moduloDivisor)
                    obj.T_reordering.running = true;
                    obj.VRux = mod(obj.VRuh,moduloDivisor);
                end
            end
            
        end
%%        
        function obj = t_ReorderingExpired(obj)
            obj.SDUsDropped = 0;
            moduloDivisor = 2^obj.sn_FieldLength;
            %update VR(UR) to the SN of the first UMD PDU with SN >= VR(UX) that has not been received
            j = 1;
            RXbuffSNs = [];
            for k = 1:length(obj.RXbuff)
                SN = readRLCPDUSN(obj,obj.RXbuff{k});
                if modulothing(SN,obj,moduloDivisor) >= modulothing(obj.VRux,obj,moduloDivisor)
                    RXbuffSNs(j,2) = SN;
                    RXbuffSNs(j,1) = modulothing(SN,obj,moduloDivisor);
                    j = j+1;
                end
            end
            RXbuffSNs = sortrows(RXbuffSNs);
            if isempty(RXbuffSNs) || RXbuffSNs(1,1) ~= modulothing(obj.VRux,obj,moduloDivisor)
                %VR(UX) not RXed 
                obj.VRur = mod(obj.VRux,moduloDivisor);
            else
                if length(RXbuffSNs(:,1)) > 1
                    for k = 1:length(RXbuffSNs(:,1))-1
                        if RXbuffSNs(k,1)+1 ~= RXbuffSNs(k+1,1)
                            obj.VRur = mod(RXbuffSNs(k,2)+1,moduloDivisor);
                            break
                        elseif k == length(RXbuffSNs) - 1
                            obj.VRur = mod(RXbuffSNs(k+1,2)+1,moduloDivisor);
                        end
                    end
                else %length(RXbuffSNs(:,1)) == 1
                    obj.VRur = mod(RXbuffSNs(1,2)+1,moduloDivisor);
                end
            end
            obj = reassemble(obj,1);
            
            obj.T_reordering.running = false;
            obj.T_reordering.timer = obj.someresetvalue;
            if modulothing(obj.VRuh,obj,moduloDivisor) > modulothing(obj.VRur,obj,moduloDivisor)
                obj.T_reordering.running = true;
                obj.T_reordering.timer = obj.someresetvalue;
                obj.VRux = mod(obj.VRuh,moduloDivisor);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % HELPER FUNCTIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function SN = readRLCPDUSN(obj,pdu)
            %Reads RLC PDU SN from header and outputs as a double
            if obj.sn_FieldLength == 5
                snBin = pdu(4:8);
            elseif obj.sn_FieldLength == 10
                snBin = pdu(7:16);
            else
                fprintf("ERROR: Bad sn_FieldLength input. sn_FieldLength Should be either 5 or 10. Simulation results will be garbage.\n")
            end
            
            SN = 0;
            pwr = length(snBin)-1;
            for k = 1:length(snBin)
                SN = SN + (2^pwr)*snBin(k);
                pwr = pwr-1;
            end
        end  
        
        function obj = reassemble(obj,mode)
            %Mode 1 (only happens when PDUs outside of window have been cleared): 
            %Reassemble RLC SDUs from any UMD PDUs with 
            %SN < updated VR(UR), remove RLC headers when doing so
            %and deliver the reassembled RLC SDUs to upper layer in 
            %ascending order of the RLC SN if not delivered before
            
            %Mode 2:
            %Reassemble RLC SDUs from any UMD PDUs with SN that falls 
            %outside of the reordering window, remove RLC headers when 
            %doing so and deliver the reassembled RLC SDUs to upper layer 
            %in ascending order of the RLC SN if not delivered before
            if ~isempty(obj.RXbuff)
                moduloDivisor = 2^obj.sn_FieldLength;
                %extract sequence numbers and buffer index from buffer, store in matrix
                RXbuffSNs = zeros(length(obj.RXbuff),3);
                for k = 1:length(obj.RXbuff)
                    %RXbuff SN
                    RXbuffSNs(k,3) = readRLCPDUSN(obj,obj.RXbuff{k});      
                    %RXbuff SN index
                    RXbuffSNs(k,2) = k;
                    %Modulothing of RXbuff SN
                    RXbuffSNs(k,1) = modulothing(RXbuffSNs(k,3),obj,moduloDivisor);
                end

                %Sort pdus by Modulothing of RXbuff SN and keep buffer index mapping too
                RXbuffSNs = sortrows(RXbuffSNs);

                %Sort RXbuff by SN
                tempBuff = cell(size(obj.RXbuff));
                for k = 1:length(obj.RXbuff)
                   tempBuff{k} = obj.RXbuff{RXbuffSNs(k,2)};
                end
                obj.RXbuff = tempBuff;

                if (mode ~= 1) && (mode ~= 2)
                    fprintf("ERROR: Select either 1 or 2 for mode input. Difference between modes is explained in the function code.\n")
                    return
                end
                l = 1; %SDU index counter
                k = 1;
                while true
                    if (mode == 1 && RXbuffSNs(k,1) < modulothing(obj.VRur,obj,moduloDivisor)) ...
                            || (mode == 2 && modulothing(obj.VRuh,obj,moduloDivisor) < RXbuffSNs(k,1))

                       %Parse SDUs, reassemble if needed
                       [FI, singlePDUSDUs] = parsePDU(obj,k);
                       %check if SDU fragment from previous reassembly is present
                       if FI(1) && obj.SDUfragmentSN >= 0 ... 
                               && modulothing(obj.SDUfragmentSN + 1,obj,moduloDivisor) == RXbuffSNs(k,1)
                           %Check if still a fragment
                           if length(singlePDUSDUs) == 1 && FI(2)
                               obj.SDUfragment = [obj.SDUfragment singlePDUSDUs{1}];
                               obj.SDUfragmentSN = RXbuffSNs(k,3);
                           else
                               %SDU is now complete, delete fragment, save whole
                               SDUs{l} = [obj.SDUfragment singlePDUSDUs{1}];
                               l = l+1;
                               obj.SDUfragment = [];
                               obj.SDUfragmentSN = -1;
                           end
                       elseif FI(1) == 0
                           %Either this is a stanadlone SDU or the start of
                           %a fragmented SDU OR starts with standalone SDU
                           %and concatenates a segmented SDU >.>
                           if FI(2) && length(singlePDUSDUs) == 1
                               %It's the start of a fragmented SDU that
                               %took up entire space of PDU
                               if ~isempty(obj.SDUfragment)
                                   obj.SDUsDropped = obj.SDUsDropped+1;
                               end
                               obj.SDUfragment = singlePDUSDUs{1};
                               obj.SDUfragmentSN = RXbuffSNs(k,3);
                           else
                               %There's no fragment to deal with
                               SDUs{l} = singlePDUSDUs{1};
                               l = l+1;
                               %delete any fragment stored in entity
                               if obj.SDUfragmentSN >= 0
                                   obj.SDUfragment = [];
                                   obj.SDUfragmentSN = -1;
                               end
                           end
                       else
                           %SDU fragments are irrecoverable; delete fragment
                           %stored in entity, skip fragment in singlePDUSDUs{1}
                           obj.SDUfragment = [];
                           obj.SDUfragmentSN = -1;
                           obj.SDUsDropped = obj.SDUsDropped+1;
                       end
                       %Now loop through singlePDUSDUs starting at 2nd index (first
                       %index was dealt with above), end at penultimate index
                       for m = 2:length(singlePDUSDUs)-1
                           SDUs{l} = singlePDUSDUs{m};
                           l = l+1;
                       end

                       %now need to check if last SDU in singlePDUSDUs is a
                       %fragment (if last SDU isn't also the first)
                       if length(singlePDUSDUs) > 1
                           if FI(2)
                               %yes this is a fragment
                               obj.SDUfragment = singlePDUSDUs{length(singlePDUSDUs)};
                               obj.SDUfragmentSN = RXbuffSNs(k,3);
                           else
                               %this is not a fragment
                               SDUs{l} = singlePDUSDUs{length(singlePDUSDUs)};
                               l = l+1;
                           end
                       end
                       %Remove PDU from buffer and row from RXbuffSNs, adjust counter
                       obj.RXbuff(k) = [];
                       RXbuffSNs(k,:) = [];
                       k = k - 1;
                    end
                    k = k + 1;
                    if k > length(obj.RXbuff)
                        break
                    end
                end
                if l == 1
                    %didn't have to reassemble anything
                    SDUs = [];
                end
                if ~isempty(SDUs)
                    for k = 1:length(SDUs)
                        obj.SDUs{length(obj.SDUs)+1} = SDUs{k};
                    end
                end
            end
        end
        
        function [FI, SDUs] = parsePDU(obj,PDUindex)
            %This function returns the FI field of a UM RLC PDU and an
            %array of all the RLC SDUs contained within it
            if obj.sn_FieldLength == 5
                %Read header this way
                
                %Get FI field
                FI = obj.RXbuff{PDUindex}(1:2);
                
                %check fixed header E field
                if obj.RXbuff{PDUindex}(3)
                    %There's an extended header
                    %Find how many SDUs are in PDU
                    numSDUs = 2;
                    Eindex = 21;
                    while obj.RXbuff{PDUindex}(Eindex)
                        %There's another SDU
                        numSDUs = numSDUs + 1;
                        Eindex = Eindex + 12;
                    end
                    %Calculate starting index of data field
                    if mod(numSDUs,2)
                        %numSDUs is odd
                        dataIndex = (2.5 + 1.5*numSDUs - 1)*8 + 1;
                    else
                        %numSDUs is even
                        dataIndex = (2 + 1.5*numSDUs - 1)*8 + 1;
                    end
                    %Parse SDUs
                    LIindex = 10;
                    LIfieldLength = 11;
                    SDUs = cell(1,numSDUs);
                    for k = 1:numSDUs
                       %extract LI field for each SDU, convert to decimal
                       LIbin = obj.RXbuff{PDUindex}(LIindex:LIindex+LIfieldLength-1);
                       LIindex = LIindex + LIfieldLength + 1;
                       LIdec = 0;
                       for j = 1:length(LIbin)
                           LIdec = LIdec + LIbin(j)*2^(length(LIbin) - j);
                       end
                       %finally ready to parse
                       SDUs{k} = obj.RXbuff{PDUindex}(dataIndex:dataIndex+LIdec-1);
                       dataIndex = dataIndex + LIdec;
                    end
                else
                    %Just one SDU after the fixed header
                    SDUs{1} = obj.RXbuff{PDUindex}(9:length(obj.RXbuff{PDUindex}));
                end
            elseif obj.sn_FieldLength == 10
                %Read header that way
                FI = obj.RXbuff{PDUindex}(4:5);
                
                %check fixed header E field
                if obj.RXbuff{PDUindex}(6)
                    %There's an extended header
                    %Find how many SDUs are in PDU
                    numSDUs = 2;
                    Eindex = 29;
                    while obj.RXbuff{PDUindex}(Eindex)
                        %There's another SDU
                        numSDUs = numSDUs + 1;
                        Eindex = Eindex + 12;
                    end
                    %Calculate starting index of data field
                    if mod(numSDUs,2)
                        %numSDUs is odd
                        dataIndex = (2.5 + 1.5*numSDUs)*8 + 1;
                    else
                        %numSDUs is even
                        dataIndex = (2 + 1.5*numSDUs)*8 + 1;
                    end
                    %Parse SDUs
                    LIindex = 18;
                    LIfieldLength = 11;
                    SDUs = cell(1,numSDUs);
                    for k = 1:numSDUs
                       %extract LI field for each SDU, convert to decimal
                       LIbin = obj.RXbuff{PDUindex}(LIindex:LIindex+LIfieldLength-1);
                       LIindex = LIindex + LIfieldLength + 1;
                       LIdec = 0;
                       for j = 1:length(LIbin)
                           LIdec = LIdec + LIbin(j)*2^(length(LIbin) - j);
                       end
                       %finally ready to parse
                       SDUs{k} = obj.RXbuff{PDUindex}(dataIndex:dataIndex+LIdec-1);
                       dataIndex = dataIndex + LIdec;
                    end
                else
                    %Just one SDU after the fixed header
                    SDUs{1} = obj.RXbuff{PDUindex}(17:length(obj.RXbuff{PDUindex}));
                end
            else
                fprintf("ERROR: bad sn_FieldLength")
            end
        end
        
        function sn = modulothing(sn,obj,moduloDivisor)
           %moduloDivisor is the Divisor you'd use when setting a state
           %variable to something, i.e. 2^sn_FieldLength
           modbase = obj.VRuh - obj.UM_Window_Size;
           sn = mod(sn-modbase,moduloDivisor);
        end
    end
end

