%This random number generator and random-variate generation functions are
%VBA translations of the C programs found in
%Law, A. M. and Kelton, W. D., ``Simulation Modeling and Analysis'',
%Singapore: The McGraw-Hill Book Co, pp. 430--431.

%Function to generate lognormal variate by transforming a normal
%Note MeanPrime and VariancePrime as the DESIRED mean and variance for the lognormal

function y = Lognormal(MeanPrime, VariancePrime, Stream)
 
    Mean = log(MeanPrime ^ 2 / sqrt(MeanPrime ^ 2 + VariancePrime));
    Variance = log(1 + VariancePrime / MeanPrime ^ 2);
 
    y = exp(Normal(Mean, Variance, Stream));

end