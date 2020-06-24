% generate 3D Euclidean distances
function [D,dDdx,dDdy,dDdz] = distSetPts3D(x)

% D is a matrix giving the distances between
% the n points specified by the nx3 vector x
% D_ij = d(x_i,x_j)^2
%
% dDdx is the anti-symmetric matrix which is the derivative of D_ij with respect to x_i
% dDdx ..."... y_i

x2 = dot(x,x,2);
D = sqrt(bsxfun(@plus,x2,x2') - 2*(x*x'));
if nargout ==1, return; end

dDdx = bsxfun(@minus,x(:,1),x(:,1)')./D;
dDdy = bsxfun(@minus,x(:,2),x(:,2)')./D;
dDdz = bsxfun(@minus,x(:,3),x(:,3)')./D;

for ii = 1:size(x,1), dDdx(ii,ii) = 0; dDdy(ii,ii) = 0; dDdz(ii,ii) = 0; end
