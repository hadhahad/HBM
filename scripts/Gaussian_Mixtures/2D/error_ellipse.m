function error_ellipse(mu, covariance)
% Function to plot an error ellipse in the 2D space
%
% input:
%   mu         - mean value vector
%   covariance - covariance matrix
%
% output:
%   none   

[eigenvec, eigenval] = eig(covariance);

% Get the index of the largest eigenvector
[max_evc_ind_c, ~] = find(eigenval == max(max(eigenval)));
max_evc = eigenvec(:, max_evc_ind_c);

% Get the largest eigenvalue
max_evl = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(max_evc_ind_c == 1)
    min_evl = max(eigenval(:,2));
    min_evc = eigenvec(:,2);
else
    min_evl = max(eigenval(:,1));
    min_evc = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(max_evc(2), max_evc(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end


% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);
phi = angle;
X0  = mu(1);
Y0  = mu(2);
a   = chisquare_val*sqrt(max_evl);
b   = chisquare_val*sqrt(min_evl);

% Ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

% Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

% Ellipse rotation
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Plot the error ellipse
fig = gcf;
plot(X0,Y0,'r*');
hold on;
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'--')
hold on;

end

