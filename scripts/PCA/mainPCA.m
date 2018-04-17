clc;
clear;

syms a x real;

%% Numerical Values:
d = 1;
r_e = 0.5;
r_a = 2;
r_x = 2;

%% Laplace:

[a_res_lap,x_res_lap, V, D] = toyPCALaplace();
a_res_lap = subs(a_res_lap);
x_res_lap = subs(x_res_lap);
V = subs(V);
V = subs(V, [a, x], [a_res_lap, x_res_lap])
D = subs(D);
D = subs(D, [a, x], [a_res_lap, x_res_lap])

toyPCALapVis(d, r_e, r_a, r_x, a_res_lap, x_res_lap);

%% Variational Bayes:
[a_res_vb,x_res_vb, s_a, s_x] = toyPCABayes(d, r_e, r_a, r_x, 0.0001);

toyPCAVBVis(d, r_e, r_a, r_x, a_res_vb, x_res_vb, s_a, s_x);