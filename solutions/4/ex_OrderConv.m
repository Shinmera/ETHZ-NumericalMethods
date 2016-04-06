clear all;
% different definitions of f with the same zero:
f = @(x) x.*exp(x)-1;   
% f = @(x) exp(x)-1./x;     % f = @(x) x-exp(-x);
x0 = 1;
x_star = fzero(f,x0);

x = x0; upd = 1;
while (abs(upd) > eps)
    fx = f(x(end));   % only 2 evaluations of f at each step
    if fx ~= 0;
        upd = fx^2 / (f(x(end)+fx)-fx);
        x = [x, x(end)-upd];
    else upd = 0;
    end
end
residual = f(x);
err      = abs(x-x_star);
log_err  = log(err);
ratios   = (log_err(3:end)-log_err(2:end-1))...
    ./(log_err(2:end-1)-log_err(1:end-2));

disp('x, err, residuals,ratios')
[x',err',residual',[0;0;ratios']] 