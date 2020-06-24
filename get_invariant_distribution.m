% compute the invariante distribution
% P is th etransition probability from i to j
function v=get_invariant_distribution(P)
	[V,D] = eig(P');
	D=diag(D);
	[~,i]=min(abs(D-1));
	assert(abs(D(i)-1) < 0.001);
	v = V(:,i);
	assert(norm(P'*v - v)/norm(v) < 1e-3);
	if v(1)<0
		v=-v;
	end
	assert(min(v)>0);
	v=v/sum(v);
end
