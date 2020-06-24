% get the Fiedler vector associated with Fan Chung's Laplacian
% a walk is given by v*P
function [Fiedler_vec,lam] = fiedler_FC(P)

	% Find the invariant measure:
	%v = pagerank(P',.95,1e-12);
	v=get_invariant_distribution(P);

	Phi = diag(v.^(1/2));
	Phiinv = diag(v.^(-1/2));
	L = eye(size(P,1))- 1/2*(Phi*P*Phiinv + Phiinv*(P')*Phi);
	max(max(abs(L - L')))/max(max(abs(L)))
	L = max(L,L');
	[V,D] = eig(L);
	lam = diag(D);
	Fiedler_vec = V(:,2);
end
