function runDevice3App(event,varargin)
    global UDs
    global baseComfortD3
    global PriceReadingD3
    global demandReductionD3
    global currentDemandD3
    global currentDemandFunk
    global gClk
    
    switch event
        case 1
            %Send demand reduction availability update
            %instead of random, do power reading
            lengthPfield = 14;
            currentDemandD3 = varargin{1};
            reductionAvailability = currentDemandD3 - demandReductionD3 - baseComfortD3;
            if reductionAvailability < 0
                reductionAvailability = 0;
            end
            ra = num2binArray(round(reductionAvailability),lengthPfield);
            currentDemand = num2binArray(round(currentDemandD3),lengthPfield);
            if isempty(UDs(3).pktNumber)
                UDs(3).pktNumber = 1;
            else
                UDs(3).pktNumber = UDs(3).pktNumber + 1;
            end
            rlcSDU = [zeros(1,77), ra, currentDemand, randi([0,1],1,2000-105)]; %
            rlcSDU = configMetaDataHeader(UDs(3),rlcSDU);

            UDs(3).RLCtxEntity.rxRLCSDU(rlcSDU);

        case 2
            %Recieved message(s)
            for k = 1:length(UDs(3).rxBuffer)
                sdu = UDs(3).rxBuffer{k};
                PriceReadingD3 = binArray2num(sdu(78:91));
                %check if price is high enough to self-regulate demand
                if PriceReadingD3 > 400
                   demandReductionD3 = currentDemandFunk(gClk/1000) - baseComfortD3;
                else
                   demandReductionD3 = 0;
                end
                rxedDemandReduction = binArray2num(sdu(106:119));
                
                if rxedDemandReduction > demandReductionD3
                    if rxedDemandReduction > currentDemandFunk(gClk/1000) - baseComfortD3
                        demandReductionD3 = currentDemandFunk(gClk/1000) - baseComfortD3;
                    else
                        demandReductionD3 = rxedDemandReduction;
                    end
                end
            end
            
        otherwise
    end
end