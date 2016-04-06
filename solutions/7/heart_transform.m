function heart_transform ()

xy = heart();
n = size(xy,2) - 1;

% rotation matrix
theta = pi/3;
A_1 = [cos(theta) sin(theta); -sin(theta) cos(theta)];
xy_1 = A_1 * xy;

% shear transformation matrix
lambda = 3;
A_2 = [1 lambda; 0 1];
xy_2 = A_2 * xy;

%evaluation points:
tt = 0:0.005:1;

% equidistant parametrization of [0,1]:
t_eq = (0:n)/n;
[pol spl pch] = curveintp (xy,t_eq,tt);
[pol_A_1 spl_A_1 pch_A_1] = curveintp (xy_1,t_eq,tt);
[pol_A_2 spl_A_2 pch_A_2] = curveintp (xy_2,t_eq,tt);

display('% === Errors for equidistant parametrization:');
fprintf('Error |pol_A_1 - A_1_pol| = %f\n', norm(pol_A_1- A_1*pol));
fprintf('Error |spl_A_1 - A_1_spl| = %f\n', norm(spl_A_1- A_1*spl));
fprintf('Error |pch_A_1 - A_1_pch| = %f\n', norm(pch_A_1- A_1*pch));
fprintf('Error |pol_A_2 - A_2_pol| = %f\n', norm(pol_A_2- A_2*pol));
fprintf('Error |spl_A_2 - A_2_spl| = %f\n', norm(spl_A_2- A_2*spl));
fprintf('Error |pch_A_2 - A_2_pch| = %f\n', norm(pch_A_2- A_2*pch));

% segment length parametrization of [0,1]:
t_seg = segment_param(xy);
[pol spl pch] = curveintp (xy,t_seg,tt);
[pol_A_1 spl_A_1 pch_A_1] = curveintp (xy_1,t_seg,tt);
t_seg = segment_param(A_2*xy); % arc-length parametrization needs to be adjusted
[pol_A_2 spl_A_2 pch_A_2] = curveintp (xy_2,t_seg,tt);

display('% === Errors for segment length parametrization:');
fprintf('Error |pol_A_1 - A_1_pol| = %f\n', norm(pol_A_1- A_1*pol));
fprintf('Error |spl_A_1 - A_1_spl| = %f\n', norm(spl_A_1- A_1*spl));
fprintf('Error |pch_A_1 - A_1_pch| = %f\n', norm(pch_A_1- A_1*pch));
fprintf('Error |pol_A_2 - A_2_pol| = %f\n', norm(pol_A_2- A_2*pol));
fprintf('Error |spl_A_2 - A_2_spl| = %f\n', norm(spl_A_2- A_2*spl));
fprintf('Error |pch_A_2 - A_2_pch| = %f\n', norm(pch_A_2- A_2*pch));

end

% segment length parametrization of [0,1]:
function t_seg = segment_param (xy)
  increments = sqrt(sum(diff(xy,1,2).^2));
  t_seg = cumsum(increments);
  t_seg = [0,t_seg/t_seg(end)];
end
