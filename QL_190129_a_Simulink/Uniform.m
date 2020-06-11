%This random number generator and random-variate generation functions are
%VBA translations of the C programs found in
%Law, A. M. and Kelton, W. D., ``Simulation Modeling and Analysis'',
%Singapore: The McGraw-Hill Book Co, pp. 430--431.
%Function to generate Uniform(Lower, Upper) variates via inverse cdf

function y=Uniform(Lower, Upper, Stream)
    y = Lower + (Upper - Lower) * lcgrand(Stream);
end