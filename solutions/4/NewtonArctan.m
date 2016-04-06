clear all; 
x0 = 2;   % initial guess
r = 1;
while (r > eps)
    x1 = x0-( 2*x0-(1+x0^2)*atan(x0) ) / (1-2*x0*atan(x0));
    % x1 = (-x0+(1-x0^2)*atan(x0))/(1-2*x0*atan(x0));
    r = abs ( (x1 - x0) / x1 );
    x0 = x1;
    fprintf ( 'x0 =  %16.14e , accuracy = %16.14e \n ' ,  x1, r );
end

figure;
x1 = x0-atan(x0)*(1+x0^2);   x2 = x1-atan(x1)*(1+x1^2);
X=[-2:0.01:2];
plot(X, atan(X),'k',...
    X,  2*(X)-(1+(X).^2).*atan((X)),'r--',...
    [x0, x1, x1, x2, x2], [atan(x0), 0, atan(x1), 0, atan(x2)],...
    [x0,x1],[0,0],'ro',[-2,2], [0,0],'k','linewidth',2); 
legend('arctan', 'g', 'Newton critical iteration');axis equal;
print -depsc2 'ex_NewtonArctan.eps'