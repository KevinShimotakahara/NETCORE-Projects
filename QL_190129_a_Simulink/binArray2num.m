function num = binArray2num(binArray)
    num = 0;
    for k = 1:length(binArray)
        num = num + binArray(length(binArray)-k+1)*2^(k-1);
    end
end