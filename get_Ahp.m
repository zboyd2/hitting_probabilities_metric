% get A^(hp,beta) from the paper
% P(i,j) is the transition probability from i to j
function Aht=get_Ahp(P,beta)

	if nargin < 2
		beta=.5;
	end

	Q = HittingTimes_L3(P); % get hitting probabilities

	% Find the invariant measure:
	%G=digraph(P);
	%v=centrality(G,'pagerank','FollowProbability',1,'Tolerance',1e-5,'MaxIterations',2000,'Importance',double(G.Edges.Weight));
	%v = pagerank(P',1,1e-2);
	n = size(P,1);
	opts.tol=1e-3;
	opts.v0=ones(n,1);
	opts.disp=0;
	opts.maxit = 100000;
	[v,~] = eigs(double(P'),1,1+1e-6,opts);
	v = abs(v) / norm(v, 1);

	% Construct the symmetric adjacency matrix M:
	if beta == .5
		Aht = spdiags(v.^(.5),0,n,n)*Q*spdiags(v.^-.5,0,n,n);
	elseif beta==1
		Aht = diag(v)*Q;
	end
	clear Q;

	assert(min(min(Aht)) >= 0);
	%max(max(abs(Aht-Aht')))/max(max(abs(Aht)));

	Aht = (Aht + Aht')/2;
	Aht = Aht - diag(diag(Aht));

end
