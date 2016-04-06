function y = arrowmatvec(d,a,x)
% Multiplying a vector with the product of two ``arrow matrices''
% Arrow matrix is specified by passing two column vectors a and d
if (length(d) ~= length(a)), error ('size mismatch'); end
% Build arrow matrix using the MATLAB function diag()
A = [diag(d(1:end-1)),a(1:end-1);(a(1:end-1))',d(end)];
y = A*A*x; 
