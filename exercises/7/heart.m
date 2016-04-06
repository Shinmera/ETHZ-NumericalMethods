function xy_heart = heart ()

% This is a fine grid of points drawing a heart:
x = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 84.5 83.5 82.5 81 79 75 70 64 60 55 50 45 40 35 30 25 20 15 10 5 0];
y = [-100.000 -94.118 -88.235 -82.353 -76.471 -70.588 -64.706 -58.824 -52.941 -47.059 -41.176 -35.294 -29.412 -23.529 -17.647 -11.765 -5.882 0 5 10 15 20 25 30 35 40 42 44 45 46 45.5 45 44 42 40 37.5 34 30 25];

% We use the symmetries:
x_heart_full = [x -x(end-1:-1:1)];
y_heart_full = [y  y(end-1:-1:1)];

% Now we choose some non-equidistant grid. You can make your own choice.
choice = [1 6 7 12 13 14 18 22 27 30 33 36 40 45 52 57 58 61 75 1];
x_heart = x_heart_full(choice);
y_heart = y_heart_full(choice);

xy_heart = [x_heart; y_heart];
