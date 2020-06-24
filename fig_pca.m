% replicate PCA figure (it depends a fair amount on random chance through the kmeans and random graph generation
N=100; % node count
nhat=3; % community count
pin=.4; % in-edge probability
pout=.1; % out-edge probability
directed=true;

addpath('util');

A = planted_partition(N,nhat,pin,pout,directed); % make the graph
n = round(N/nhat); % community size

if ~is_connected(A)
	msgID = 'dhp_kmedoids:NotConnected';
	msg = 'Graph is not connected';
	e = MException(msgID,msg);
	throw(e);
end

% get distances
P = diag(sum(A,2))^(-1)*A; % random walk is v*P not P*v
global dhp;
dhp = -log(get_Ahp(P));
for i=1:N; dhp(i,i)=0; end

% pca
c=pca(dhp);

% ground truth partition (for coloring)
gt = ones(N,1);
for a=2:nhat
	gt( (1+(a-1)*n) : (a*n) ) = a;
end
gt(end)=nhat;

% plot
scatter(c(:,1),c(:,2),[],gt,'filled')
axis square off
