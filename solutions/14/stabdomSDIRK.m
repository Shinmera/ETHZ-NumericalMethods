function stabdomSDIRK(gamma)
% STABDOMSDIRK Plot the stability domain for the two-stage SDIRK Method
%   (Singularly Diagonal Implicit Runge-Kutta). GAMMA is the parameter on the
%   diagonal of the matrix A.
%
%   Butcher Table:
%
%      gamma   |   gamma        0
%    1 - gamma | 1 - 2*gamma  gamma
%    _______________________________
%              |    1/2        1/2
%   
% For Example:
% stabdomSDIRK(1-sqrt(2)/2);   % 2nd order, l-stable.
% stabdomSDIRK((3-sqrt(3))/6); % 3rd order, not a-stable, not l-stable.
% stabdomSDIRK((3+sqrt(3))/6); % 3rd order, a-stable, not l-stable.
%
% See also SDIRK, SDIRKSTEP.

if nargin < 1
gamma = (3+sqrt(3))/6; % default gamma
end

% generate grid
[X, Y] = meshgrid(-15:0.1:15);
Z = X + i*Y;

% define stability function
S = (1 + Z*(1 - 2*gamma) + Z.^2*(0.5 - 2*gamma + gamma^2)) ./ (1-gamma*Z).^2;

% plot stability domain
figure
contourf(X,Y,abs(S),0:0.1:1);
colormap(hot);
colorbar;

title(['Stability Domain for 2-Stage SDIRK, \gamma = ' num2str(gamma) ])
xlabel('Re')
ylabel('Im')
axis square
grid on