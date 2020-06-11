%This random number generator and random-variate generation functions are
%VBA translations of the C programs found in
%Law, A. M. and Kelton, W. D., ``Simulation Modeling and Analysis'',
%Singapore: The McGraw-Hill Book Co, pp. 430--431.
%Function to generate a random integer


function y = Random_integer(prob_distrib, Stream)

 % Generate a U(0,1) random variate.
    U = lcgrand(Stream);
    
%  Return a random integer in accordance with the (cumulative) distribution
%  function prob_distrib.
    y = 1;
    while (U >= prob_distrib(y))
        y = y + 1;
    end
end