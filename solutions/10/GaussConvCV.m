function GaussConvCV(f_hd)
if nargin<1;    f_hd = @(t) sinh(t);    end;
g = @(t) t.*f_hd(sin(t)).*cos(t);    
I_exact = quad(@(t) asin(t).*f_hd(t),-1,1,eps)
%I_exact = quad(@(t) g(t),-pi/2,pi/2,eps)

n_max = 50;  nn = 1:n_max;
err = zeros(size(nn));
for j = 1:n_max
    [x,w] = gaussquad(nn(j));
    I = pi/2*dot(w, g(x*pi/2));
    % I = GaussArcSinCV(f_hd,nn(j));   % using pcode
    err(j) = abs(I - I_exact);
end

close all; figure; 
semilogy(nn,err,[1,n_max],[eps,eps],'--','linewidth',2);    
title('{\bf Convergence of Gauss quadrature for rephrased problem}','fontsize',14);
xlabel('{\bf n = # of evaluation points}','fontsize',14);
ylabel('{\bf error}','fontsize',14);
legend('Gauss quad.','eps');
print -depsc2 '../PICTURES/GaussConvCV.eps';
