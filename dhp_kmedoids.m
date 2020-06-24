% Auxilliary script used by fig_communities.m
function p = dhp_kmedoids(varargin)

	% parse input
	p = inputParser;
	addOptional(p,'pin',.4);
	addOptional(p,'pout',.1);
	addParameter(p,'use_pca',true);
	addParameter(p,'N',300);
	addParameter(p,'directed',true);
	addParameter(p,'nhat',3);
	parse(p,varargin{:});
	pin=p.Results.pin;
	pout=p.Results.pout;
	N=p.Results.N;
	nhat=p.Results.nhat;
	p = p.Results;

	% generate graph
	A = planted_partition(N,nhat,pin,pout,p.directed);
	n = round(N/nhat);

	% ground truth
	gt = ones(N,1);
	for a=2:nhat
		gt( (1+(a-1)*n) : (a*n) ) = a;
	end
	gt(end)=nhat;

	% confirm connectedness
	if ~is_connected(A)
		msgID = 'dhp_kmedoids:NotConnected';
		msg = 'Graph is not connected';
		e = MException(msgID,msg);
		throw(e);
	end

	% pca on A
	c=pca(A);
	ind = kmeans(c(:,1:(nhat-1)),nhat);
	%ind = kmeans(A-mean(mean(A)),nhat);
	p=zeros(3,1);
	p(1) = purityMeas(ind,gt);

	% get hitting probability distances
	P = diag(sum(A,2))^(-1)*A; % random walk is v*P not P*v
	global dhp;
	dhp = -log(get_Ahp(P));
	for i=1:N; dhp(i,i)=0; end

	% pca on hitting probability distances
	c=pca(dhp);
	ind=kmeans(c(:,1:(nhat-1)),nhat);
	p(2) = purityMeas(ind,gt);

	% use kmedoids
	ind=kmedoids((1:N)',nhat,'Distance',@dist);
	p(3) = purityMeas(ind,gt);

end

function d=dist(ZI,ZJ)
	global dhp;
	d=dhp(ZJ,ZI);
end
