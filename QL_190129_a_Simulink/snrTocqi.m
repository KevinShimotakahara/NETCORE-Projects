%% Map snr to cqi
function [cqi] = snrTocqi(snr)
    cqi_table = 1:15;
    snr_table = 0:2:28;

    [~, col] = find(snr_table <= snr);

    if(isempty(col))
        cqi = 1;
    else
        cqi = cqi_table(col(end));
    end
    
    if(cqi > 9)
        cqi = 9;
    end
end



