function [KF,ll]=kalman(KF,y, yp)
% function evaluates the Kalman filter given by KF:
%
% KF.A   - state transition matrix
% KF.Q   - state noice covariance
% KF.C   - observation matrix
% KF.R   - observation noise matrix
% KF.x   - state estimate
% KF.P   - state variance


xnn1 = KF.A*KF.x;
KF.P  = KF.A*KF.P*KF.A' + KF.Q;

%Data update
Ry = KF.C*KF.P*KF.C' + KF.R;
iRy = inv(Ry);
K = KF.P*KF.C'*iRy;
KF.P = KF.P- K*KF.C*KF.P; % P = P -KCP;
% yp=KF.C*xnn1;
KF.x = xnn1 + K*(y-yp);

if nargin>1
    ll = -0.5*log(det(Ry))- 0.5*(y-yp)'*iRy*(y-yp);
end

