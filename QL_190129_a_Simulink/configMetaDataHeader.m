function sdu = configMetaDataHeader(UD,sduIn) 
    global gClk
    sduIn(1,1:53) = [num2binArray(UD.nodeID,5), num2binArray(UD.pktNumber,24), num2binArray(gClk,24)];
    sdu = sduIn;
end

function binArray = num2binArray(num,arrayLength)
    %Used for populating binary headers with decimal values and adding any
    %padding if needed
    binArray = zeros(1,arrayLength);
    numStr = dec2bin(num);
    for k = 1:length(numStr)
        digit = str2double(numStr(length(numStr)-k+1));
        binArray(length(binArray)-k+1) = digit;
    end
end