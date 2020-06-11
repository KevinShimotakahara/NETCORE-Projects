
lambda = 50; % arrival rate per second 
T=500*1e-3; % simulation time in second (10 hours)
delta=1e-3; % simulation step size in second
N=T/delta; % number of simulation steps
event=zeros(N,1); % array recording at each step if a "packet" arrived. 
 % initialize it to zeros
R=rand(size(event)); % generate a random array (with elements in [0,1]) of the same size as "event"
event(R<lambda*delta)=1; % set each element of event to 1 with probability lambda*delta
inds=find(event==1); % getting indices of arrivial
int_times=diff(inds)*delta; % interarrival times in seconds
edges=0:1e-3:0.5; % define histogram bin
count=histc(int_times,edges);

figure; bar(edges,count,'histc'); % draw histogram of absolute count