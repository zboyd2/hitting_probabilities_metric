% computes generalized effective resistance
function R = extEffRes(A)

% inputs:
% A is the adjacency matrix for a (possibly) directed graph
% A_{i,j} = 1 if (i,j) is an edge and = 0 otherwise
%
% outputs: 
% R is the extended effective resistance as defined in 
% Young, Scardovi, and Leonard (2013)
% note: R^{1/2} is a metric, but R is not
% 
% The notation used below is from this paper

n = size(A,1);
d_out = sum(A,2); % out degrees
L = diag(d_out) - A; % satisfies  L*ones(n,1) = 0

% Pi and Q satisfy Q*Pi = Q
Pi = eye(n) - (1/n)*ones(n); % projection matrix
[V,~] = eig(Pi); 
Q = V(:,2:end)';

% reduced Laplacian
rL = Q*L*Q';

% solve the Lyapunov equation for Sig
Sig = lyap(rL,-eye(n-1)); 

% X plays the role of L^{-1}
X = 2*Q'*Sig*Q;

R = zeros(n);
for ii = 1:n
    for jj = ii:n
        r = X(ii,ii) + X(jj,jj) - 2*X(ii,jj);        
        R(ii,jj) = r;
        R(jj,ii) = r;
    end
end

