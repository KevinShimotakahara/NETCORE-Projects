function y = advanceLTEsim(u)
    %put Medhat's main loop here; eventually will have outputs that will
    %influence smart grid as well
    %y will be a vector of event outputs depending on the SDUs RXed from
    %various D2D entities during this TTI   
    global gClk

    Main()

    global UDs
    global ackTimer2
    global ackTimer3
    global ackTimer4
    global demandReductionD2
    global demandReductionD3
    global demandReductionD4
%     try
    %Increment timeout timers (if active)
    if ~isempty(ackTimer2) && ackTimer2 >= 0
        ackTimer2 = ackTimer2 + 1;
    end
    if ~isempty(ackTimer3) && ackTimer3 >= 0
        ackTimer3 = ackTimer3 + 1;
    end
    if ~isempty(ackTimer4) && ackTimer4 >= 0
        ackTimer4 = ackTimer4 + 1;
    end
    
    %call app1 if ack timeout
    if (~isempty(ackTimer2) && ackTimer2 > 700) || (~isempty(ackTimer3) && ackTimer3 > 700) || (~isempty(ackTimer4) && ackTimer4 > 700)
        runDevice1App(2)
    end
        
    y = zeros(10,1);
    if ~isempty(UDs(1).rxBuffer)
        if isempty(UDs(1).rxBuffer{1})
            UDs(1).rxBuffer(1) = [];
            if ~isempty(UDs(1).rxBuffer)
                runDevice1App(1);
                UDs(1).rxBuffer = [];
                UDs(1).rxBuffer = cell(1);
            else
                UDs(1).rxBuffer = cell(1);
            end
        else
            runDevice1App(1);
            UDs(1).rxBuffer = [];
            UDs(1).rxBuffer = cell(1);
        end
    end
    if ~isempty(UDs(2).rxBuffer)
        if isempty(UDs(2).rxBuffer{1})
            UDs(2).rxBuffer(1) = [];
            if ~isempty(UDs(2).rxBuffer)
                runDevice2App(2)
                UDs(2).rxBuffer = [];
                UDs(2).rxBuffer = cell(1);
            else
                UDs(2).rxBuffer = cell(1);
            end
        else
            runDevice2App(2)
            UDs(2).rxBuffer = [];
            UDs(2).rxBuffer = cell(1);
        end
    end
     if ~isempty(UDs(3).rxBuffer)
        if isempty(UDs(3).rxBuffer{1})
            UDs(3).rxBuffer(1) = [];
            if ~isempty(UDs(3).rxBuffer)
                runDevice3App(2)
                UDs(3).rxBuffer = [];
                UDs(3).rxBuffer = cell(1);
            else
                UDs(3).rxBuffer = cell(1);
            end
        else
            runDevice3App(2)
            UDs(3).rxBuffer = [];
            UDs(3).rxBuffer = cell(1);
        end
    end
     if ~isempty(UDs(4).rxBuffer)
        if isempty(UDs(4).rxBuffer{1})
            UDs(4).rxBuffer(1) = [];
            if ~isempty(UDs(4).rxBuffer)
                runDevice4App(2)
                UDs(4).rxBuffer = [];
                UDs(4).rxBuffer = cell(1);
            else
                UDs(4).rxBuffer = cell(1);
            end
        else
            runDevice4App(2)
            UDs(4).rxBuffer = [];
            UDs(4).rxBuffer = cell(1);
        end
     end
     
     %update DRevent Signals
     y(8) = demandReductionD2;
     y(9) = demandReductionD3;
     y(10) = demandReductionD4;
%     catch
%         bleh = 1;
%     end
end
