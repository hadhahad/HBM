clc;
clear all; %#ok<CLALL>

% Initial parameters:
mu1 = [2 1];
mu2 = [-2 2];
mu3 = [0 -1];
Sig1 = eye(2);
Sig2 = eye(2);
Sig3 = diag([2, 0.1]);
N = [300,300,400];

plot_truepdf = false;
plot_sample = true;

%% Plot True Probability Density
if plot_truepdf
    x = -2:0.1:2;
    y = -2:0.1:2;
    [X_axis,Y_axis] = meshgrid(x, y);
    F1 = mvnpdf([X_axis(:) Y_axis(:)], mu1, Sig1);
    F1 = reshape(F1,length(Y_axis),length(X_axis));
    F2 = mvnpdf([X_axis(:) Y_axis(:)], mu2, Sig2);
    F2 = reshape(F2,length(Y_axis),length(X_axis));
    F3 = mvnpdf([X_axis(:) Y_axis(:)], mu3, Sig3);
    F3 = reshape(F3,length(Y_axis),length(X_axis));
    figure(1);
    surf(x,y,N(1)/sum(N).*F1);
    hold on;
    surf(x,y,N(2)/sum(N).*F2);
    hold on;
    surf(x,y,N(3)/sum(N).*F3);
    xlabel('x'); ylabel('y'); zlabel('Probability Density');
end

%% Random Samples
rng default  % For reproducibility
R = [mvnrnd(mu1,Sig1,N(1)); mvnrnd(mu2,Sig2,N(2)); mvnrnd(mu3,Sig3,N(3))];
if plot_sample
    figure(2);
    plot(R(:,1),R(:,2),'+');
end

%% Inference
K = 3;
[mu,sig] = mixture_gauss2D(R,K,4);

%% Visualization
figure(2);
hold on;
plot(mu(:,1),mu(:,2),'r*')