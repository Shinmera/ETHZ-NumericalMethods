%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Some visualuzation (no subproblem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Visualize the complex plane, the path, the ellipse and the set S
thetas = linspace(0,2*pi);
rho = 1.5;
max_rho = exp(asinh(pi/3))

gamma = @(rho, theta) cos(theta - 1i*log(rho));
f = @(z) 1 ./ (1 + exp(-3*z));

n1 = 3;
S = 2/3*pi*(-(n1):(n1-1)) + pi / 3;
fS = exp(-3*S*1i) + 1; % Check if S*i is actually a zero of exp(-3z) +1

I = [-1,1];

hold on
plot(I,[0,0],'r.-','LineWidth',2);

Gamma = gamma(rho,thetas);
plot(real(Gamma), imag(Gamma),'g-');
Gamma = gamma(max_rho,thetas);
plot(real(Gamma), imag(Gamma),'g--');

plot(0,S,'b.');

legend('[0,1]','\gamma_\delta, \delta = 1.5','maximal \gamma_\delta, \delta = 2.4952', 'S')
title('Visualization of the complex plane')
xlabel('real part')
ylabel('imagnary part')
print -depsc '../PICTURES/ellipse.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subproblem (iv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 1000; % Discretization parameter
[m, rho] = cheby_approx(max_rho, N, f, gamma) % Approximate M

degrees = 1:20; % Test for degrees 1 to 20

% Theoretical upper bound computed for (iv)
%upper_bound = 2..^(-degrees)./(rho.^(degrees+1)-1)*m*sqrt(2*(1./max_rho^2 + max_rho^2)) / 2;
upper_bound = 2../(rho.^(degrees+1)-1)*m*sqrt(2*(1./max_rho^2 + max_rho^2)) / 2;

% Compute the actual error for any degree
errInf = [];
for n = degrees

    % Cheby nodes and value of f at the nodes
    k = 0:n;
    t = cos((2*k + 1) / (2*(n+1)) * pi);
    y = f(t);
    
    % Evaluate the error at a couple of points
    x = linspace(-1,1);

    F = f(x);
    LnF = intpolyval(t,y,x);
    %LnF = polyval(polyfit(t,y,n),x);
    
    
    errInf = [errInf, max(abs(F - LnF))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Some visualization (no subproblem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Visualize Chebychev nodes
figure
plot(t,0,'bo');
title('Chebychev nodes')
xlabel('t')
print -depsc '../PICTURES/chebnodes.eps'

% Visualize f and interpolation of f, with deg 20
figure
plot(x,F,'r.--', x, LnF,'b-');
title('Function and its interpolant')
xlabel('t')
legend('f(t)', 'L_{20}(t)');
print -depsc '../PICTURES/interpolant.eps'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subproblem (v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subproblem (v): comparison of our upper ound and acutal error
figure
semilogy(degrees, errInf, degrees, upper_bound);
title('Real and upper bound for L^\infty error')
xlabel('degree')
ylabel('L^\infty error')
legend('Computed error','Upper bound')
print -depsc '../PICTURES/error.eps'
