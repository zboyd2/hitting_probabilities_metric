% compute the hitting PROBABILITIES (I know the function name is bad) of the row-stochastic matrix M. Slow (O(n^4)) version.
function [Q] = HittingTimes(F)
% Input a stochastic matrix, output the matrix of probs. of leaving i 
% then hitting j before returning to i
% O(N^4) time

N = size(F,1);
Q = zeros(N,N);
A = zeros(N,N);
Aj = zeros(N,N);
A = eye(N) - F;

for j = 1:N

    Aj = A;
    ej = zeros(1,N);
    ej(j) = 1;
    Aj(j,:) = ej;

    Ajinv = inv(Aj);

    for i = 1:N
	if i ~= j
        	Q(i,j) = Ajinv(i,j)/Ajinv(i,i);
	end
    end
end






    % Increase Accuracy with Jacobi Iterations:
%     Dj = diag(diag(Aj));
%     Rj = Aj-Dj;
%     errit = 1;
%     while errit > 1e-12
%         Ajinv0 = Dj\[eye(N) - Rj*Ajinv];
%         errit = max(max(abs(Ajinv - Ajinv0)));
%         Ajinv = Ajinv0;
%     end
