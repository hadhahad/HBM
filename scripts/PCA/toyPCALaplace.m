function [a_res, x_res, V, D] = toyPCALaplace()

syms d a x real;
syms r_e r_x r_a positive;

log_pax_d = -(d-a*x)^2/(2*r_e) - a^2/(2*r_a) - x^2/(2*r_x);
Gradient = gradient(log_pax_d,[a,x]);
[a_hat,x_hat] = solve(Gradient,[a,x]);
a_hat = simplify(a_hat);
x_hat = simplify(x_hat);
Hessian = -simplify(hessian(log_pax_d, [a, x]));
sigma_res = inv(Hessian)
[V, D] = eig(sigma_res);

% non-complex positive values
a_res = eval(a_hat(2,1));     
x_res = eval(x_hat(2,1));

end