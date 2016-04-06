function [pol spl pch] = curveintp (xy,t,tt)

x = xy(1,:);
y = xy(2,:);

% compute the slopes at the extremes (necessary for complete spline):
x_start = (x(2)-x(1)) ./ (t(2)-t(1));
y_start = (y(2)-y(1)) ./ (t(2)-t(1));
x_end = (x(end)-x(end-1)) ./ (t(end)-t(end-1));
y_end = (y(end)-y(end-1)) ./ (t(end)-t(end-1));
    
% compute the interpolants
polx = intpolyval(t,x,tt);
poly = intpolyval(t,y,tt);
  
% complete splines, using the extended vector
splx = spline(t,[x_start x x_end],tt);
sply = spline(t,[y_start y y_end],tt);
pchx = pchip(t,x,tt);
pchy = pchip(t,y,tt);

pol = [polx; poly];
spl = [splx; sply];
pch = [pchx; pchy];

end
