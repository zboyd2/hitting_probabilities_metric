% Get fiedler vector and eigenvalue of the hitting probability distance using P
% P(i,j) is the transition probability from i to j
function [Fiedler_vec,lam] = ht_spectral(P)
% Build Hitting Time Matrix:

Aht = get_Ahp(P);

% Get Fiedler vector and eigenvalue
L = diag(sum(Aht)) - Aht;
[V,lam] = eigs(double(L),2,-1e-6);
lam=diag(lam);
Fiedler_vec = V(:,2)';
end
