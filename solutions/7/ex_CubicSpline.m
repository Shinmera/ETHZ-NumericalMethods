function ex_CubicSpline(alpha,beta)
if nargin==0;  alpha=-1; beta=-11/3;  end;

x1 = -1:0.01:0;
y1 = (x1+1).^4+alpha*(x1-1).^4+1;
x2 = 0:0.01:1;
y2 = -x2.^3-8*alpha*x2+1;
x3 = 1:0.01:2;
y3 = beta*x3.^3 + 8.*x3.^2+11/3;

x = [x1 x2 x3]; y = [y1 y2 y3];
nodes = [-1 0 1 2];
data = [y1(1) y1(end) y2(end) y3(end)];

close all;
plot(x,y,nodes,data,'ro','linewidth',2);
legend('cubic spline', 'data points','Location','SouthEast');
xlabel('x','fontsize',14); ylabel('s(x)','fontsize',14);
title('Cubic spline with parameters','fontsize',14)
print -depsc2 'ex_CubicSpline.eps'