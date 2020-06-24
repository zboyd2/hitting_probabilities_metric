% Get fiedler vector and eigenvalue of a naive symmetrization of P
% P(i,j) is the transition probability from i to j
function [Fiedler_vec,lam] = naive_fiedler(P)
	%A = (P+ P')/2;
	A=max(P,P');
	L = diag(sum(A)) - A;
	[V,D] = eig(L);
	lam = diag(D);
	Fiedler_vec = V(:,2);
end
