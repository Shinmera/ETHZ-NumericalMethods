C               = 2;
p               = 1.5;
eps_max         = C^(1/(1-p));

ngp             = 100; % number of grid points
eps_lin         = linspace(0, eps_max, ngp);
tau_lin         = linspace(0, eps_max, ngp);
[eps_msh,tau_msh] = meshgrid(eps_lin(2:(end-1)),tau_lin(2:(end-1)));

kmin = @(eps, C, p, tau) ceil(log( ( log(tau) + (1/(p-1)).*log(C) ) ./ log(C^(1/(p-1)) .* eps)  ) ./ log(p) );
k = kmin(eps_msh, C, p, tau_msh);

% Consider only gridpoints where: eps larger as tau
for ne = 1:ngp-2
    for nt = 1:ngp-2
        if (ne > nt)
           k(ne,nt) = 0;
        end
    end
end

% Plotting
pcolor(eps_msh,tau_msh,k)
colorbar()
title('Minimal number of iterations for error < \tau')
xlabel('\epsilon_0')
ylabel('\tau')
xlim([0,eps_max])
ylim([0,eps_max])
shading flat