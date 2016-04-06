function y = arrowmatvec2(d,a,x)
if (length(d) ~= length(a)), error ('size mismatch'); end
% Recursive matrix-vector multiplication to obtain A*A*x = A*(A*x)
y = A(d,a,A(d,a,x));
end

% Efficient multiplication of a vector with the ``arrow matrix''
function Ax = A(d,a,x)
Ax = d.*x;
Ax(1:end-1) = Ax(1:end-1) + a(1:end-1)*x(end);
Ax(end) = Ax(end) + a(1:end-1)'*x(1:end-1);
end
