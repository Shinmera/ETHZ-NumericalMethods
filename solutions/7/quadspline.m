% part d of quadratic splines problem
% given: t = nodes (n-1 vect.),
%        y = data (n vect.),
%        x = evaluation pts (N vect.)
% create the interpolating quadratic spline and evaluate in x
function eval_x = quadspline(t,y,x)
  % the number of nodes:
  n = length(y);         % N = length(x);
  % ensure nodes and data are line vectors:
  t = t(:)'; y = y(:)';
  % create (n+3) extended vectors using the periodicity:
  ext_t = [t(end)-1,0,t,1,1+t(1)];
  % increments in t:
  de_t = diff(ext_t);                      % (n+2)
  dde_t = ext_t(3:end) - ext_t(1:end-2);   % (n+1)
  
  % build the three n-vectors that define the matrix
  A = de_t(2:n+1) ./ (2*dde_t(1:n));
  B = 1 + de_t(3:n+2) ./ (2*dde_t(2:n+1)) + ...
      de_t(1:n) ./ (2*dde_t(1:n));
  C = de_t(2:n+1) ./ (2*dde_t(2:n+1));
  % Assembly of the matrix, be careful with spdiags and transpose!
  M = spdiags([C; B; A].', (-1:1), n,n).';
  M(1,n) = A(1);
  M(n,1) = C(n);
  
  % solve the system for the coefficients c:
  c = (M\(y(:))).';

  % compute the coefficients d_1,...,d_n:
  ext_c = [c,c(1)];
  d = 2*( ext_c(2:end).*de_t(2:n+1)+ c.*de_t(3:n+2) )...
      ./dde_t(2:n+1);
  ext_d = [d(n),d];
  % evaluate the interpolating spline in x:
  eval_x = zeros(1,length(x));
  
  % loop over the n intervals
  ind_end = 0;
  for i=1:n
      left_t = ext_t(i+1);
      right_t = ext_t(i+2);
      % find the points in x in the interval [t(i-1), t(i)):
      ind_start = ind_end+1;
      if x(ind_start) >= left_t
          while ind_end < length(x) && x(ind_end+1) < right_t
              ind_end = ind_end + 1;
          end
      end
      ind_eval = ind_start:ind_end;
      tau = (x(ind_eval) - left_t)./( right_t-left_t );
      eval_x(ind_eval) = d(i)*tau.^2 + ...
          c(i)*4*tau.*(1-tau) + ext_d(i)*(1-tau).^2;
  end
end