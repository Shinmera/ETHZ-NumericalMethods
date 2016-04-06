function heart_Sol ()

xy = heart();
n = size(xy,2) - 1;

%evaluation points:
tt = 0:0.005:1;

figure;
hold on;

% equidistant parametrization of [0,1]:
t_eq = (0:n)/n;
[pol spl pch] = curveintp (xy,t_eq,tt);
subplot(1,2,1); xlabel('Equidistant Parametrization');
plot_interpolations (xy,pol,spl,pch);

% segment length parametrization of [0,1]:
t_seg = segment_param(xy);
[pol spl pch] = curveintp (xy,t_seg,tt);
subplot(1,2,2); xlabel('Segment Length Parametrization');
plot_interpolations (xy,pol,spl,pch);

hold off;
print -depsc2 '../PICTURES/ex_CurveIntp.eps'

end

% segment length parametrization of [0,1]:
function t_seg = segment_param (xy)
  increments = sqrt(sum(diff(xy,1,2).^2));
  t_seg = cumsum(increments);
  t_seg = [0,t_seg/t_seg(end)];
end

% plotting function
function plot_interpolations (xy,pol,spl,pch)
  plot(xy(1,:),xy(2,:),'o', pol(1,:),pol(2,:),'-.', ...
       spl(1,:),spl(2,:),'-', pch(1,:),pch(2,:),'--', 'linewidth',2);
  axis equal;
  axis([-105 105 -105 105]);
  legend('data','polynomial','spline','pchip','Location','Southoutside');
end