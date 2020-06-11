%%

function metadata = readMetadata(sdu,varargin)
    if isempty(varargin)
        sdu = cell2mat(sdu);
    end
    metadata.nodeID = binArray2num(sdu(1:5));
    metadata.packetNum = binArray2num(sdu(6:29));
    metadata.TOC = binArray2num(sdu(30:53));
    metadata.TOTX = binArray2num(sdu(54:77));
    
    function num = binArray2num(binArray)
        num = 0;
        for k = 1:length(binArray)
            num = num + binArray(length(binArray)-k+1)*2^(k-1);
        end
    end
end


