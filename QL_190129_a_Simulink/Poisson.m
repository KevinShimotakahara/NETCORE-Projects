
% PLS: Payload Size (in bits)

function [TxFlag] = Poisson(Prob, lambda)
        
    x = randi(10, 1);
    pd = makedist('Poisson', lambda);
    ngen = pdf(pd, x);
    
    if(ngen > Prob)
        TxFlag = 1;
    else
        TxFlag = 0;
    end
end
