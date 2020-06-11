
function CI = ConfInt(x)

SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.05  0.95],length(x)-1);        % T-Score
CI = mean(x) + ts*SEM;                      % Confidence Intervals

end


