function  y_ev = PWlineIntp (t,y,t_ev)
% compute and evaluate piecewise linear interpolant
% t and y     data vector of the same size
% t_ev        vector with evaluation points
% --> y_ev    column vector with evaluations in t_ev

t=t(:);  y = y(:);  t_ev = t_ev(:);     
n = size(t,1)-1;    % # intervals
if n~=size(y,1)-1; error('t, y must have the same size'); end;

[t,I] = sort(t);    % sort t and y if not sorted
y = y(I);

y_ev = zeros(size(t_ev));
for k=1:n
    t_left  = t(k);
    t_right = t(k+1);
    ind = find ( t_ev >= t_left & t_ev < t_right );
    y_ev(ind) = y(k) + (y(k+1)-y(k))/(t_right-t_left)*(t_ev(ind)-t_left);
end
% important! take care of last node:
y_ev(find(t_ev == t(end))) = y(end);