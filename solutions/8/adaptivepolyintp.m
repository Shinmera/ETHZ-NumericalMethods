function [t,errs] = adaptivepolyintp(f,a,b,tol,N)
% Adaptive polynomial interpolation of function f:[a,b] -> R.
% argument 'f' must be a handle that allows vectorized evaluation,
% argument 'tol' specifies the relative tolerance,
% the optional argument 'N' passes the number of sampling points.
if (nargin < 5), N = 1000; end
sp = (a:(b-a)/N:b);      % N+1 equdistant sampling points
fvals = f(sp);           % evaluate function at sampling points
fmax  = max(abs(fvals)); % Maximum of f in sampling points
t = sp(floor(N/2));      % set of interpolation nodes
y = fvals(floor(N/2));   % function values at interpolation nods
errs = [];
for i=1:N
  err = abs(fvals - intpolyval(t,y,sp)); % modulus of interpolation error 
  [mx,idx] = max(err); errs = [errs, mx]; 
  % finishes, once maximal pointwise interpolation error below threshold
  if (mx < tol*fmax), return; end
  t = [t,sp(idx)];
  y = [y,fvals(idx)];
end
error('Desired accuracy could not be reached');
