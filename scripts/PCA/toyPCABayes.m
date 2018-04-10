function [a_res, x_res, s_a, s_x] = toyPCABayes(d, r_e, r_a, r_x, eps)

% Initial conditions:
a_hat = 6;
x_hat = 2.5;
var_a = 0.25;
var_x = 0.25;

% Iterations:
for i = 1:1:10000
    a_sec_mom = a_hat^2 + var_a;
    x_sec_mom = x_hat^2 + var_x;
    
    var_a = ((x_sec_mom / r_e) + (1 / r_a))^(-1);
    var_x = ((a_sec_mom / r_e) + (1 / r_x))^(-1);
    
    aux_a = x_hat * d * var_a / r_e;
    aux_x = a_hat * d * var_x / r_e;
    if (abs(aux_a - a_hat) < eps && abs(aux_x == x_hat) < eps)
        break;
    else
        a_hat = aux_a;
        x_hat = aux_x;
    end 
end

a_res = a_hat;
x_res = x_hat;
s_a = sqrt(var_a);
s_x = sqrt(var_x);

end