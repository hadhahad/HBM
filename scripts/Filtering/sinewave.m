freq = 50*2*pi;
dt = 1/1000;
% gen data
N = 200;
vt = (1:N)*dt*freq;
if 0
    % step change
    a_true = 5
    y0 = [sin(freq*dt*(1:100)) a_true.*sin(freq*dt*(1:100)+1.5)];
else
    % slowly varying
    a_true = (1+1*(1:100)/100);
    y0 = [sin(freq*dt*(1:100)) a_true.*sin(freq*dt*(1:100))];
end
y = y0 + 0.2*randn(size(y0));

res=zeros(2,200);
 KF.A=[1,0;0,1];
 KF.Q=10^(-1)*[1,0;0,1];
 KF.R=0.1;
 KF.x=[0;0.1];
 KF.P=eye(2);
for i=1:1:200,
    
    KF.C=[sin(res(2,i) + freq*dt*i), i*cos(res(2,i) + freq*dt*i)];
    res(:,i)=KF.x;
    y_est(1,i)=res(1,i)*sin(freq*dt*i+res(2,i));
    
    [KF,ll]=kalman(KF,y(1,i));

end

fig=figure(1)
subplot(2,1,1)
hold off
plot(y0)
hold on
plot(y,'.');

% TODO estimate a and phi by EKF

