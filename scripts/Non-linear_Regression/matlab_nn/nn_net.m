clear all
syms w11 w12 w13 w14 w15 w16
syms b11 b12 b13 b14 b15 b16
syms w21 w22 w23 w24 w25 w26
syms b2
syms y x

% first layer
z = [tanh(w11*x+b11);...
    tanh(w12*x+b12);...
    tanh(w13*x+b13);...
    tanh(w14*x+b14);...
    tanh(w15*x+b15);...
    tanh(w16*x+b16);];

yp= w21*z(1) + w22*z(2) + w23*z(3) + w24*z(4) + w25*z(5) + w26*z(6) + b2;

loss = (y-yp)^2;

allw = [w11 w12 w13 w14 w15 w16 b11 b12 b13 b14 b15 b16 w21 w22 w23 w24 w25 w26 b2];
g = gradient(loss,allw.')

z = matlabFunction(z);
gr = matlabFunction(g);
nnz = @(x,W,b)z(b(1),b(2),b(3),b(4),b(5),b(6),W(1),W(2),W(3),W(4),W(5),W(6),x);
nng = @(x,y,W,b,W2,b2)gr(b2,b(1),b(2),b(3),b(4),b(5),b(6),...
    W(1),W(2),W(3),W(4),W(5),W(6),...
    W2(1),W2(2),W2(3),W2(4),W2(5),W2(6),...
    x,y);
nny = @(z,W2,b2) (W2(1)*z(1) + W2(2)*z(2) + W2(3)*z(3) + W2(4)*z(4) + W2(5)*z(5) + W2(6)*z(6) + b2);
nnyp = @(x,W1,b1,W2,b2) nny(nnz(x,W1,b1),W2,b2);
save nnz nny nng nnyp