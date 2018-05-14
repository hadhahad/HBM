% train
load nn_data
%P = [W1 b1 W2 b2];
P = 0.1*rand(3*6+1,1);
for it=1:50000
    G = zeros(3*6+1,1);
    L = 0;
    
    % compute gradient for all points
    for i = 1:length(x0)
        G = G + nng(x0(i),y(i),P(1:6),P(7:12),P(13:18),P(19));
        L = L + (y(i) - nnyp(x0(i),P(1:6),P(7:12),P(13:18),P(19)))^2;
    end
    % gradient descent
    P = P - 0.001*G;
    
    LL(it) = L;
end

%%
fig=figure(3);
hold off
semilogy(LL);
title('Convergence of the loss function')

fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0. 4 3];
fig.PaperSize = [4 3];
print('nn_L','-dpdf')


% validate        
for j = 1:length(x0)
    yp(j) = nnyp(x0(j),P(1:6),P(7:12),P(13:18),P(19));
    zp(:,j) = nnz(x0(j),P(1:6),P(7:12));
end

%% 
fig=figure(1);
hold off
plot(x0,y,'rx');
hold on
plot(x0,yp);
legend('data','NN')

fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0. 4 3];
fig.PaperSize = [4 3];
print('nn_fit','-dpdf')


%%
fig=figure(2);
hold off
plot(x0,zp);
title('basis functions')

fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0. 4 3];
fig.PaperSize = [4 3];
print('nn_z','-dpdf')
