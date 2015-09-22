function y = arrowmatvec2(d, a, x)

if (length(d) /= length(a))
    error('size mismatch');
end

A = [diag(d(1:end-1)), a(1:end-1);
     (a(1:end-1))', d(end)];

y = (A*(A*x));

end