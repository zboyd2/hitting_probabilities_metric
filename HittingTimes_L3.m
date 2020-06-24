% compute the hitting PROBABILITIES (I know the function name is bad) of the row-stochastic matrix M. Fast (O(n^3)) version.
function [Q] = HittingTimes_L3(M)
	% Input a stochastic matrix, output the matrix of probs. of leaving i 
	% then hitting j before returning to i
	% O(N^3) time

	N = size(M,1);

	prec = 'double';

	e1 = zeros(1,N,prec);
	e1(1) = 1;
	A1inv = eye(N)-M; % Not A1inv yet, but this saves on memory
	%A = A1;
	A1inv(1,:) = e1;
	A1inv = inv(A1inv);

	Q = zeros(N,N,prec);
	Q(:,1) = A1inv(:,1)./diag(A1inv);

	M = M*A1inv;

	detCj = (1 + diag(M))*(1-M(1,1)) + M(:,1).*M(1,:)';
	CjInv = zeros(2,2,N,prec);
	CjInv(1,1,:)=(1-M(1,1)) ./ detCj;
	CjInv(1,2,:)=M(:,1) ./ detCj;
	CjInv(2,1,:)=-M(1,:)' ./ detCj;
	CjInv(2,2,:)=(1+diag(M)) ./ detCj;

	M1 = zeros(N,2,N,prec);
	M1(:,1,:) = A1inv;
	M1(:,2,:) = repmat(-A1inv(:,1),1,N); % so each column is the same

	M2 = zeros(2,N,N,prec);
	M2(1,:,:) = M';
	M2(2,:,:) = repmat(M(1,:),N,1)';

	%Aj = zeros(N,N);
	Ac = zeros(N,1);
	Ad = zeros(N,1);
	for j = 2:N

		%        Aj = A;
		%        ej = zeros(1,N);
		%        ej(j) = 1;
		%        Aj(j,:) = ej;
		%        Cj = eye(2) + [ej; e1]*F*A1inv*[ej', -e1'];
		%        Ajinv = A1inv - A1inv*[ej', -e1']*(inv(Cj))*[ej; e1]*F*A1inv;
		%
		%	Q(:,j) = Ajinv(:,j)./diag(Ajinv);

		%        Ajinv = A1inv - [A1inv(:,j) -A1inv(:,1)]*CjInv(:,:,j)*[M(j,:); M(1,:)];
		%	Q(:,j) = Ajinv(:,j)./diag(Ajinv);

		assert(all(all([A1inv(:,j) -A1inv(:,1)] == M1(:,:,j))));
		assert(all(all([M(j,:); M(1,:)] == M2(:,:,j))));
		Ac = A1inv(:,j) - M1(:,:,j)*CjInv(:,:,j)*M2(:,j,j);
		Ad = diag(A1inv) - sum(M1(:,:,j)'.*(CjInv(:,:,j)*M2(:,:,j)))';
		Q(:,j) = Ac./Ad;
	end
	clear Ajinv;
	clear A1inv;
	clear M;
	clear M1;
	clear M2;
	Q = Q-diag(diag(Q));

end

