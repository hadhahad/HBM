clc;
clear;

freq = 50*2*pi;
N = 500;
N_half = N ./ 2;
dt = 1 / (N * 5);

%% Generate Data
a_mode = 3;
if a_mode == 1
    % step change
    a_true = 5;
    a_true_abs_arr(1, 1:N_half) = 1;
    a_true_abs_arr(1, (N_half+1):N) = a_true;
    phi_true_arr(1, 1:N_half) = 0;
    phi_true_arr(1, (N_half+1):N) = 1.5;
    y0 = [sin(freq*dt*(1:N_half)) a_true.*sin(freq*dt*(1:N_half)+1.5)];
elseif a_mode == 2
    % slowly varying
    a_true = 1 + 1*(1:N_half)/N_half;
    a_true_abs_arr(1, 1:N_half) = 1;
    a_true_abs_arr(1, (N_half+1):N) = a_true;
    phi_true_arr(1, 1:N_half) = 0;
    phi_true_arr(1, (N_half+1):N) = 0;
    y0 = [sin(freq*dt*(1:N_half)) a_true.*sin(freq*dt*(1:N_half))];
elseif a_mode == 3
    a_true = 5.*sin(1*(1:N_half)/N_half);
    a_true_abs_arr(1, 1:N_half) = 1;
    a_true_abs_arr(1, (N_half+1):N) = a_true;
    phi_true_arr(1, 1:N_half) = 0;
    phi_true_arr(1, (N_half+1):N) = pi / 2;
    y0 = [sin(freq*dt*(1:N_half)) a_true.*sin(freq*dt*(1:N_half)+pi/2)];
end
y = y0 + 0.2*randn(size(y0));

%% Perform EKF

% Parameters initialization:
curr_state = zeros(2, N);
y_hat = zeros(1,N);
KF.A = eye(2);
KF.Q = 0.01*eye(2);
KF.R = 0.01;
KF.x = [0;0];
KF.P = eye(2);

% EKF estimation
for t = 1:1:N
    curr_state(:,t) = KF.x;
    KF.C = [sin(freq * dt * t + curr_state(2,t)), ...
         curr_state(1,t) * cos(freq * dt * t + curr_state(2,t))];
    y_hat(1,t) = curr_state(1,t) * sin(freq * dt * t + curr_state(2,t));
    [KF,ll] = kalman(KF,y(1,t),y_hat(1,t));
end


%% Visualize
fig = figure(1);

subplot(2,2,[1,2]);
plot(y0);
hold on
plot(y,'.');
hold on;
plot(y_hat);
grid on
title('Sine Wave Simulation - EKF');
xlabel('Time, t');
ylabel('y');
legend('true', 'noisy','EKF');

subplot(2,2,3);
plot(a_true_abs_arr);
hold on;
plot(abs(curr_state(1,:)));
grid on;
legend('true', 'EKF');
xlabel('Time, t');
ylabel('A_{t}');

subplot(2,2,4);
plot(phi_true_arr);
hold on;
plot(abs(curr_state(2,:)));
grid on;
legend('true', 'EKF');
xlabel('Time, t');
ylabel('\phi_{t}');