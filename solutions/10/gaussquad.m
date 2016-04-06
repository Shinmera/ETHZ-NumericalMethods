function [x,w]=gaussquad(n)
% Computation of weights and nodes of \Blue{$n$}-point Gaussian quadrature rule
% on the interval \Blue{$[-1,1]$}.
if (n==1), x = 0; w = 2;
else
  b = zeros(n-1,1); % \label{gw:1}
  for i=1:(n-1), b(i)=i/sqrt(4*i*i-1); end %\label{gw:2}
  J=diag(b,-1)+diag(b,1); [ev,ew]=eig(J); %\label{gw:3}
  x=diag(ew); w=(2*(ev(1,:).*ev(1,:)))'; %\label{gw:4}
end