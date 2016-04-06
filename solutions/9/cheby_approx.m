function [m, ret_rho] = cheby_approx(max_rho, N, f, gamma )
% Returns an approximation of min_\rho max_abs(f(z)) / d(\gamma,[-1,1])
% and the rho for which the minimum is attaiined
% max_rho: maximal value of rho, compute min in (1,max_rho)
% N: number of intervals to discretize (0,pi) and (1,max_rho)
% f: function handle for f
% gamma: function handle for gamma(rho, theta)

% Upper bound for min and return value
m = Inf;

% Discretize the values of \rho
rhos = linspace(1,max_rho,N);

% Remove 1 and exp(asinh(pi/3)) (interval is open)
rhos = rhos(2:end-1);
ret_rho = rhos(1);

% Discretize the ellipse
thetas = linspace(0, 2*pi, N);

    
% Loop over a bunch of delta
for rho = rhos,
    % z discretize the values of \gamma_\delta
    z = gamma(rho,thetas);
    
    % Sup of f on \gamma
    M = max(abs(f(z)));
    
    % Distance is the one between -1 and the minimal real part of the ellipse (i.e. semi-major-axis - 1):
    d = (rho + 1/rho)/2 - 1;
    
    % If you do not belive on the semi-major-axis distance formula:
    %d = min( abs(imag(z)).*(abs(real(z))<=1) + abs(z-1).*(real(z)>1) + abs(z+1).*(real(z)<-1));
    
    if M / d < m,
        ret_rho = rho;
    end
    m = min(m, M / d);
end

end