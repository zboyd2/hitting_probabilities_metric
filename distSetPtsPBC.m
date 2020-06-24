% generate 2D distances with periodic boundary conditions
function [Dc,dDcdx,dDcdy] = distSetPtsPBC(x,W,H)

% D is a matrix giving the distances between
% the n points specified by the nx2 vector x
% D_ij = d(x_i,x_j)^2
%
% distances are calculated on the Tori (periodic bc)
%
% dDdx is the anti-symmetric matrix which is the derivative of D_ij with respect to x_i
% dDdx ..."... y_i
%
% W and H are width and height of the domain

n = size(x,1);
%L=5; % note: domain is [-L,L]
%W = 10; H = 10*sqrt(3)/2; %width and height of domain
%W = 10; H = 10; %width and height of domain

Dxy = @(x,y) sqrt(bsxfun(@plus,dot(x,x,2),dot(y,y,2)') - 2*(x*y'));

D = zeros(n,n,9);
D(:,:,1) = Dxy(x,x);
D(:,:,2) = Dxy(x+[zeros(n,1),H*ones(n,1)],x);
D(:,:,3) = Dxy(x+[zeros(n,1),-H*ones(n,1)],x);
D(:,:,4) = Dxy(x+[W*ones(n,1),zeros(n,1)],x);
D(:,:,5) = Dxy(x+[-W*ones(n,1),zeros(n,1)],x);
D(:,:,6) = Dxy(x+[W*ones(n,1),H*ones(n,1)],x);
D(:,:,7) = Dxy(x+[-W*ones(n,1),-H*ones(n,1)],x);
D(:,:,8) = Dxy(x+[W*ones(n,1),-H*ones(n,1)],x);
D(:,:,9) = Dxy(x+[-W*ones(n,1),H*ones(n,1)],x);

[d,i] = min(reshape(D,n^2,9),[],2);
Dc = reshape(d,n,n); I = reshape(i,n,n);
% D contains the tori distances,
% I keeps track of how the distance was measured

if nargout ==1, return; end

Dx=zeros(n,n); Dy=zeros(n,n); % x and y differences 
for ii = 1:9
    if ii==1, px = zeros(n,1); py = zeros(n,1); end
    if ii==2, px = zeros(n,1); py = H*ones(n,1); end
    if ii==3, px = zeros(n,1); py = -H*ones(n,1); end
    if ii==4, px = W*ones(n,1); py = zeros(n,1); end
    if ii==5, px = -W*ones(n,1); py = zeros(n,1); end
    if ii==6, px = W*ones(n,1); py = H*ones(n,1); end
    if ii==7, px = -W*ones(n,1); py = -H*ones(n,1); end 
    if ii==8, px = W*ones(n,1); py = -H*ones(n,1); end
    if ii==9, px = -W*ones(n,1); py = H*ones(n,1); end

    temp = bsxfun(@minus,x(:,1)+px,x(:,1)');
    Dx(I==ii) = temp(I==ii);
    temp = bsxfun(@minus,x(:,2)+py,x(:,2)');
    Dy(I==ii) = temp(I==ii);
end

dDcdx = Dx./Dc; dDcdy = Dy./Dc;
for ii = 1:n, dDcdx(ii,ii) = 0; dDcdy(ii,ii) = 0; end
