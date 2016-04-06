function plot_adaptivepolyintp ()
f1 = @(t) sin(exp(2*t));
f2 = @(t) sqrt(t) ./ ( 1 + 16*t.^2 );
a = 0; b = 1; tol = 1e-6;
[~,err1] = adaptivepolyintp (f1,a,b,tol);
[~,err2] = adaptivepolyintp (f2,a,b,tol);

semilogy(1:length(err1), err1, 'r-', 1:length(err2), err2, 'b-');
legend('f_1', 'f_2');
title('{\bf Error decay for adaptive polynomial interpolation}');
xlabel('{\bf number of nodes}');
ylabel('{\bf error}');

print -depsc '../PICTURES/plot_adaptivepolyintp.eps'

