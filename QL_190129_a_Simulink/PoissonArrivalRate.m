clear
clc

Seed = 2^5;
rng(Seed);

TTI     = 1e-3;
lambda	= 1/(2e-3);     % arrival rate per second (1 packet per 2 milliseconds)
T		= 500 * 1e-3;   % simulation time in second (500 milliseconds)
delta	= 0.1e-3; 		% simulation step size in second
N		= T/delta; 		% number of simulation steps
M       = TTI/delta;    % number of simulation steps per TTI
event	= zeros(N,1); 	% array recording at each step if a "packet" arrived. 

 % initialize it to zeros
R			= rand(size(event)); 	% generate a random array (with elements in [0,1]) of the same size as "event"
event(R<lambda*delta) = 1; 			% set each element of event to 1 with probability lambda*delta
inds		= find(event==1); 		% getting indices of arrivial
int_times	= diff(inds)*delta*1e3; % interarrival times in milli-seconds

event2 = reshape(event, M, []);
Nevent = sum(event2, 1);            % number of packets arriving each TTI
times       = inds*delta*1e3;       % arrival times in milli-seconds









