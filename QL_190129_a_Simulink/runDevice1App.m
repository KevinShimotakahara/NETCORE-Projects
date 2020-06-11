function runDevice1App(event,varargin)
    %This is the smart grid device's application layer. Based on the event
    %experienced, it will make a decision, which usually ends up in a
    %communication message being sent to the RLC layer for processing.
    global UDs
    global ackTimer2
    global ackTimer3
    global ackTimer4
    global powerReadingD2
    global powerReadingD3
    global powerReadingD4
    global gClk
    global p
    global reductionAvailabilityD2
    global reductionAvailabilityD3
    global reductionAvailabilityD4
    global reductions
    
    switch event
        case 1
            %Recieved a message
            for k = 1:length(UDs(1).rxBuffer)
                sdu = UDs(1).rxBuffer{k};
                metadata = readMetadata(sdu,1);
                if length(UDs(1).rxBuffer{k}) == 77
                    %Is an Ack
                    switch metadata.nodeID
                        case 2
                            ackTimer2 = -1;
                        case 3
                            ackTimer3 = -1;
                        case 4
                            ackTimer4 = -1;
                        otherwise
                    end
                else
                    %Is Power Report
                    P = binArray2num(UDs(1).rxBuffer{k}(92:105));
                    ra = binArray2num(UDs(1).rxBuffer{k}(78:91));
                    switch metadata.nodeID
                        case 2
                            powerReadingD2 = P;
                            reductionAvailabilityD2 = ra;
                        case 3
                            powerReadingD3 = P;
                            reductionAvailabilityD3 = ra;
                        case 4
                            powerReadingD4 = P;
                            reductionAvailabilityD4 = ra;
                        otherwise
                    end
                end
            end
        case 2
            %ACK timeout
            DRevent()
        case 3
            %Routine price info messaging and SS power reading

            powerSecondary = varargin{1};
            if powerSecondary < -5000
                %Demand response event initiated
                reductionTarget = abs(5000+powerSecondary)+300;
                %reductions = [0,0,0];
                availabilities = [reductionAvailabilityD2,reductionAvailabilityD3,reductionAvailabilityD4;1,2,3];
                while ~isempty(availabilities)
                    [minAvailability,index] = min(availabilities(1,:));
                    if ceil(reductionTarget/length(availabilities(1,:))) < floor(minAvailability)
                        for k = 1:length(availabilities(1,:))
                            reductions(availabilities(2,k)) = reductions(availabilities(2,k)) + ceil(reductionTarget/length(availabilities(1,:)));
                        end
                        break
                    else
                        for k = 1:length(availabilities(1,:))
                            reductions(availabilities(2,k)) = reductions(availabilities(2,k)) + floor(minAvailability);
                            reductionTarget = reductionTarget - floor(minAvailability);
                            availabilities(1,k) = availabilities(1,k) - floor(minAvailability);
                        end
                        availabilities(:,index) = [];
                    end
                end
            elseif norm(reductions) > 0
                %must continue to manage demands of houses
                sumReductions = 0;
                for k = 1:length(reductions)
                    sumReductions = sumReductions + reductions(k);
                end
                if powerSecondary - sumReductions > -5000
                    %can leave DR mode
                    reductions = [0,0,0];
                else
                    %cannot leave DR mode, but can reduce reductions
                    reductionTarget = abs(4700+(powerSecondary - sumReductions));
                    splitdatshit = ceil(reductionTarget/3);
                    reductions = [splitdatshit,splitdatshit,splitdatshit];
                    %reductions = [0,0,0];
%                     availabilities = [reductionAvailabilityD2,reductionAvailabilityD3,reductionAvailabilityD4;1,2,3];
%                     while ~isempty(availabilities)
%                         [minAvailability,index] = min(availabilities(1,:));
%                         if ceil(reductionTarget/length(availabilities(1,:))) < floor(minAvailability)
%                             for k = 1:length(availabilities(1,:))
%                                 reductions(availabilities(2,k)) = ceil(reductionTarget/length(availabilities(1,:)));
%                             end
%                             break
%                         else
%                             for k = 1:length(availabilities(1,:))
%                                 reductions(availabilities(2,k)) = reductions(availabilities(2,k)) + floor(minAvailability);
%                                 reductionTarget = reductionTarget - floor(minAvailability);
%                                 availabilities(1,k) = availabilities(1,k) - floor(minAvailability);
%                             end
%                             availabilities(:,index) = [];
%                         end
%                     end
                end
            end
            
            lengthFields = 14;
            Price = num2binArray(round(p(gClk)),lengthFields);
            reductionD2 = num2binArray(round(reductions(1)),lengthFields);
            reductionD3 = num2binArray(round(reductions(2)),lengthFields);
            reductionD4 = num2binArray(round(reductions(3)),lengthFields);
            
            if isempty(UDs(1).pktNumber)
                UDs(1).pktNumber = 1;
            else
                UDs(1).pktNumber = UDs(1).pktNumber + 1;
            end

            rlcSDU = [zeros(1,77), Price, reductionD2, reductionD3, reductionD4, randi([0,1],1,2000-77-4*lengthFields)]; %
            rlcSDU = configMetaDataHeader(UDs(1),rlcSDU);

            UDs(1).RLCtxEntity.rxRLCSDU(rlcSDU);

        otherwise
    end
end
function DRevent()
    global ackTimer2
    global ackTimer3
    global ackTimer4
    global UDs
    
    if isempty(ackTimer2) || ackTimer2 >= 0
        ackTimer2 = 0;
    end
    if isempty(ackTimer3) || ackTimer3 >= 0
        ackTimer3 = 0;
    end
    if isempty(ackTimer4) || ackTimer4 >= 0
        ackTimer4 = 0;
    end
    
    if isempty(UDs(1).pktNumber)
        UDs(1).pktNumber = 1;
    else
        UDs(1).pktNumber = UDs(1).pktNumber + 1;
    end
    
    rlcSDU = [zeros(1,77), randi([0,1],1, 28*9 - 77)];
    rlcSDU = configMetaDataHeader(UDs(1),rlcSDU);

    UDs(1).RLCtxEntity.rxRLCSDU(rlcSDU);

end