%% Cancellation
h = 10.^(-20:0);
x = 1.2;

% Derivative
g1 = (sin(x+h) - sin(x)) ./ h; % Naive
g2 = 2 .* cos(x + h * 0.5) .* sin(h * 0.5) ./ h; % Better
ex = cos(x); % Exact

% Plot
loglog(h, abs(g1-ex), 'r',h, abs(g2-ex), 'b', h, h, 'k--');
title('Error of the approximation of f''(x_0)');
legend('g_1','g_2', 'O(h)');
xlabel('h');
ylabel('| f''(x_0) - g_i(x_0,h) |');
grid on
