function toyPCAVBVis(d, r_e, r_a, r_x, a_res, x_res, sigma_a, sigma_x)

a_line = -2.5:0.05:2.5;
x_line = -2.5:0.05:2.5;
t = 0:0.1:2 * pi + 0.1;

pr_density = zeros(100);
for i = 1:1:100
    for j = 1:1:100
        a = a_line(i);
        x = x_line(j);
%         pr_density(i,j) = sqrt(2*pi*r_e)^(-1) * exp(-(d-a*x)^2/(2*r_e));
        pr_density(i,j) = exp(-(d-a*x)^2/(2*r_e) - ...
            a^2/(2*r_a) - x^2/(2*r_x));
    end
end


figure;

X = meshgrid(a_line(1:end-1), x_line(1:end-1));
contour(X, X', pr_density);
hold on;
plot(a_res + sigma_a * cos(t), x_res + sigma_x * sin(t), 'm');
hold on;
plot(a_res, x_res, '*k');

xlabel('a');
ylabel('x');
title('Toy PCA - VB, d=1');
legend('p(d|a,x)', 'solution variability', 'solution');
grid on;
% colorbar;

end