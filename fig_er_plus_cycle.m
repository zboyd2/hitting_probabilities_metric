% replicate er graph plus cycle from paper
N=20; % size of er graph without the cycle
c_length=8; % cycle length
%N=714; % size of er graph without the cycle
%c_length=286; % cycle length
p=0.5; % connection density of the er graph
prec='single'; % low precision saves memory if you want to go big

addpath('util')

% A_ij is the directed edge weight from i to j
A = zeros(N+c_length,prec);
A(1:N,1:N) = double(rand(N) < p); % er part

% cycle
B=zeros(c_length,prec);
for i=1:c_length-1
	B(i,i+1)=1;
end
B(c_length,1)=1;
A(N+1:end,N+1:end) = B;
clear B;

% bidirectional connection
A(N,N+1)=1;
A(N+1,N)=1;

% Add in connections with the rest of the graph (one direction only)
for i=1:c_length
	for l=1:2*(round(N*p)-1)
		A(randi(N),N+i)=3;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(A,1);
P = diag(sum(A,2))^(-1)*A; % random walk is v*P not P*v
clear A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get fiedler vectors and eigenvalues
[Fiedler_vec, lam] = naive_fiedler(P);
[Fiedler_vecFC, lamFC] = fiedler_FC(P);
[Fiedler_vecHT,lamHT] = ht_spectral(P); % hitting probabilities
