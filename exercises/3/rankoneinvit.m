function lmin = rankoneinvit(d,tol)
if (nargin < 2), tol = 1E-6; end
ev = d; 
lmin = 0.0;
lnew = min(abs(d));

while (abs(lnew-lmin) > tol*lmin)
   lmin = lnew;
   M = diag(d) + ev*ev';
   ev = M\ev;
   ev = ev/norm(ev);
   lnew = ev'*M*ev;
end
lmin = lnew;
