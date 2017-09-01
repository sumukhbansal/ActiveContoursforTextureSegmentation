function [  out ] = log_pd(p,q)
% Computes log map on manifold of positive definite matrices.
[U,Dp]=eig(double(p));
g=U*(Dp)^(.5);
y=(g^-1)*q*(g^-1)';
[V,Dq]=eig(double(y));
Dq=ones(size(Dq))+ Dq-eye(size(Dq));    % intermediate step
out=(g*V)*log(Dq)*(g*V)'; % how to find log map of diagonal matrix
end

