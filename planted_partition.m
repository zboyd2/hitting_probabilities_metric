% generate (possibly directed) planted partition graph
% N -- node count
% nhat -- community count
% pin/out -- in/out probabilities
% directed (bool) -- chould the graph be directed?
function A=planted_partition(N,nhat,pin,pout,directed)
	n = round(N/nhat);
	A = rand(N);
	for a=1:nhat
		for b=1:nhat
			inda = (1+(a-1)*n) : (a*n);
			indb = (1+(b-1)*n) : (b*n);
			if a==nhat
				inda = (1+(a-1)*n) : N;
			end
			if b==nhat
				indb = (1+(b-1)*n) : N;
			end

			if a==b
				p = pin;
			else
				p = pout;
			end
			A(inda,indb)  = A(inda,indb) < p;
		end
	end
	if ~directed
		A = A - tril(A,-1) + triu(A,1)';
	end
	A=A-diag(diag(A));
end
