clc;
clear all; %#ok<CLALL>

% Initial parameters:
mu1 = [1 1];
mu2 = [-2 2];
mu3 = [0 -1];
Sig1 = eye(2);
Sig2 = eye(2);
Sig3 = diag([2, 0.1]);
N = [300,300,400];

plot_truepdf = true;
plot_solution = true;

%% Plot True Probability Density
if plot_truepdf
    x = -4:0.1:4;
    y = -4:0.1:4;
    [X_axis,Y_axis] = meshgrid(x, y);
    F1 = mvnpdf([X_axis(:) Y_axis(:)], mu1, Sig1);
    F1 = reshape(F1,length(Y_axis),length(X_axis));
    F2 = mvnpdf([X_axis(:) Y_axis(:)], mu2, Sig2);
    F2 = reshape(F2,length(Y_axis),length(X_axis));
    F3 = mvnpdf([X_axis(:) Y_axis(:)], mu3, Sig3);
    F3 = reshape(F3,length(Y_axis),length(X_axis));
    figure(2);
    contour(x,y,N(1)/sum(N).*F1);
    hold on;
    contour(x,y,N(2)/sum(N).*F2);
    hold on;
    contour(x,y,N(3)/sum(N).*F3);
    grid on;
    xlabel('x'); ylabel('y'); zlabel('Probability Density');
end

%% Random Samples
rng default  % For reproducibility
R = [mvnrnd(mu1,Sig1,N(1)); mvnrnd(mu2,Sig2,N(2)); mvnrnd(mu3,Sig3,N(3))];

%% Inference
K = 3;
[mu,sig,al] = mixture_gauss2D(R,K,4);
str = sprintf('%f ', al); 
fprintf('Cluster membership probabilities: %s\n', str);

%% Visualization
if plot_solution
    figure(1);
    hold on;
    plot(R(:,1),R(:,2),'.');
    for i = 1:K
        error_ellipse(mu(i,:),sig(:,:,i))
    end
    xlabel('x'); ylabel('y');
    grid on;
    title('EM - Gaussian Mixture in 2D');
    legend('simulation', 'cluster mean','variance');
end