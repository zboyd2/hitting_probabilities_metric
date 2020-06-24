% auxilliary script to support fig_distance_compare.m. It prepares gd, which is a collection of graphs/graph generators
gd = {pt_torus,pt_torus_hole,pt_H,pt_circle,pt_sphere,@square_lattice};
function [D,A]=square_lattice()
	N = 10;
	dx = 1/N;
	x = [0:dx:1-dx]';
	N1 = length(x)^2;
	e = ones(N1,1);
	
	[X,Y] = meshgrid(x);
	
	X = X(:);
	Y = Y(:);
	
	xsl = [X,Y];
	
	[D,dDdx,dDdy] = distSetPtsPBC(xsl,1,1);
	D=real(D);

	%% Square Lattice:
	N = 10;
	e = ones(N,1);
	Id = spdiags(e, 0, N, N);
	A1 = spdiags([e e e e], [-N+1, -1, 1, N-1], N, N);
	A = kron(Id,A1) + kron(A1,Id);
end

%% Point Cloud on the flat torus:
function D=pt_torus()
n=25^2; % number points
W=2*pi; H=2*pi; % domain is [0,W]x[0,H] volume is W*H = n 
x0 = rand(n,2); x0(:,1)=W*x0(:,1); x0(:,2)=H*x0(:,2); x0=x0(:);

x = [x0(1:end/2),x0(end/2 + 1:end)];
x(:,1) = mod(x(:,1),W); x(:,2) = mod(x(:,2),H);

% careful: distSetPts calculates distances, not squared distances
[D,~,~] = distSetPtsPBC(x,W,H);
end

%% Point Cloud on Flat Torus with a hole:
function D=pt_torus_hole()
n=36^2; % number points
W=2*pi; H=2*pi; % domain is [0,W]x[0,H] volume is W*H = n 
x0 = 0*rand(n,2); 

% Put the hole at W/2,H/2, radius R
R = pi/2;
index = 0;
while index < n
	z = rand(1,2);
	z(1) = W*z(1);
	z(2) = H*z(2);
	if ((z(1) - W/2)^2 + (z(2)-H/2)^2) >= R
		index = index + 1;
		x0(index,1) = z(1);
		x0(index,2) = z(2);
	end
end


x0=x0(:);

x = [x0(1:end/2),x0(end/2 + 1:end)]; 
x(:,1) = mod(x(:,1),W); x(:,2) = mod(x(:,2),H);

x(1,1)
x(1,2)

% careful: distSetPts calculates distances, not squared distances
[D,~,~] = distSetPtsPBC(x,W,H);
end

%% Point Cloud on H shaped domain:
function D=pt_H()
n=36^2; % number points
W=2*pi; H=2*pi; % domain is [0,W]x[0,H] volume is W*H = n 
x0 = 0*rand(n,2); 

% Put the hole at W/2,H/2, radius R
R1 = pi/2;
R2 = pi/4;
index = 0;
while index < n
	z = rand(1,2);
	z(1) = W*z(1);
	z(2) = H*z(2);
	if abs(z(1)-W/2) >= R1 || abs(z(2)-H/2) <= R2
		index = index + 1;
		x0(index,1) = z(1);
		x0(index,2) = z(2);
	end
end
%plot(x0(:,1),x0(:,2),'o')


x0=x0(:);

x = [x0(1:end/2),x0(end/2 + 1:end)]; 
x(:,1) = mod(x(:,1),W); x(:,2) = mod(x(:,2),H);

xcoord1 = x(1,1)
xcoord2 = x(1,2)

% careful: distSetPts calculates distances, not squared distances
[D,~,~] = distSetPts(x);  % Use non-periodic distances
end

%% Point Cloud a circle:
function D=pt_circle()
n=1000; % number points
W=2*pi; H=0; % domain is [0,W]x[0,H] volume is W*H = n 

x0 = rand(n,1); x0(:,1)=W*x0(:,1); x0(:,2)=0*x0(:,1); x0=x0(:);

x0=x0(:);

x = [x0(1:end/2),x0(end/2 + 1:end)]; 
x(:,1) = mod(x(:,1),W); x(:,2) = 0*mod(x(:,2),H);

% careful: distSetPts calculates distances, not squared distances
[D,~,~] = distSetPtsPBC(x,W,10*W);
end

%% Point Cloud a sphere:
function D = pt_sphere()
n=1000; % number points
x0 = -1+2*rand(n,3);
R = sqrt(x0(:,1).^2 + x0(:,2).^2 + x0(:,3).^2); 
x0(:,1)=x0(:,1)./R; x0(:,2)=x0(:,2)./R; x0(:,3)=x0(:,3)./R;
x = [x0(:,1),x0(:,2),x0(:,3)];
[D,~,~] = distSetPts3D(x); % careful: distSetPts calculates distances, not squared distances
end

