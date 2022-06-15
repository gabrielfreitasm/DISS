clear
close all
clc

warning off

syms x1 x2 u real
vars = [x1; x2; u];

x = [x1; x2];

% system
% example 1: Valmorbida
f = [1.5*x1^2 + x2^2 - 0.5*x1^3; -x2 + x1*x2];
g = [-1; -0.5];

% example 2: Valmorbida
% f = [2*x1^3 + x1^2*x2 - 6*x1*x2^2 + 5*x2^3; 0];
% g = [0; 1];

% example: Moez
% f = [-x1+x2+x1^2+x1*x2-x1^3+x1^2*x2-x1*x2^2+2*x2^3;
%      -x1+1.5*x2-x1^2-0.5*x1*x2-x1^3-x1^2*x2+0.5*x1*x2^2-2*x2^3];
% g = [0; 1];

monh = monomials([x1; x2], [1 2 3]);
h = ones(1, length(monh))*monh;

y = h;

% initial condition
x0 = [0; 0];

prog = sosprogram(vars);

[prog, V] = sospolyvar(prog, monomials(x, [2 3 4]));
[prog, Q] = sospolymatrixvar(prog, monomials(vars, 0), [length(h), length(h)], 'symmetric');
[prog, S] = sospolymatrixvar(prog, monomials(vars, 0), [length(h), length(u)]);
[prog, R] = sospolymatrixvar(prog, monomials(vars, 0), [length(u), length(u)], 'symmetric');
% [prog, T] = sospolyvar(prog, monomials(x, 2));

gradV = [diff(V, x1); diff(V, x2)];
% hessV = [diff(gradV(1), x1), diff(gradV(1), x2); diff(gradV(2), x1), diff(gradV(2), x2)];

% A = gradV'*f - h'*Q*h + T;
A = gradV'*f - h'*Q*h;
B = 0.5*gradV'*g - h'*S;
C = -R;

M = [A, B;
    B', C];

prog = sosineq(prog, V - 10^-5*(x1^4 + x2^4));
% prog = sosineq(prog, x'*hessV*x - 10^-5*(x1^2 + x2^2));
% prog = sosineq(prog, T - 10^-5*(x1^2 + x2^2));
prog = sosmatrixineq(prog, -M, 'quadraticMineq'); 
% or alternatively using Schur complement in M

sol = sossolve(prog);

V = sosgetsol(sol, V);
[Qv, Zv] = findsos(V);

M = sosgetsol(sol, M);
[Qm, Zm, Hm] = findsos(-M);

Q = double(sosgetsol(sol, Q));
S = double(sosgetsol(sol, S));
R = double(sosgetsol(sol, R));

return

% OBS: give K to plot the phase portrait!
% phase portrait
[x1, x2] = meshgrid(-1:.1:1, -1:.1:1);

% example 1: Valmorbida
u = K*(x1.^3 + (x1.^2).*x2 + x1.^2 + x1.*(x2.^2) + x1.*x2 + x1 + x2.^3 + x2.^2 + x2);
% u = K*(x1.^4 + (x1.^3).*x2 + x1.^3 + (x1.^2).*(x2.^2) + (x1.^2).*x2 + x1.^2 + x1.*(x2.^3) + x1.*(x2.^2) + x1.*x2 + x2.^4 + x2.^3 + x2.^2);
dx1 = - 0.5*x1.^3 + 1.5*x1.^2 + x2.^2 - u;
dx2 = x1.*x2 - x2 - 0.5*u;

% example 2: Valmorbida
% u = K*(x1.^3 + (x1.^2).*x2 + x1.^2 + x1.*(x2.^2) + x1.*x2 + x1 + x2.^3 + x2.^2 + x2);
% % u = K*(x1.^4 + (x1.^3).*x2 + x1.^3 + (x1.^2).*(x2.^2) + (x1.^2).*x2 + x1.^2 + x1.*(x2.^3) + x1.*(x2.^2) + x1.*x2 + x2.^4 + x2.^3 + x2.^2);
% dx1 = 2*x1.^3 + (x1.^2).*x2 - 6*x1.*(x2.^2) + 5*x2.^3;
% dx2 = u;

% example: Moez
% u = K*(x1.^3 + (x1.^2).*x2 + x1.^2 + x1.*(x2.^2) + x1.*x2 + x1 + x2.^3 + x2.^2 + x2);
% % u = K*(x1.^4 + (x1.^3).*x2 + x1.^3 + (x1.^2).*(x2.^2) + (x1.^2).*x2 + x1.^2 + x1.*(x2.^3) + x1.*(x2.^2) + x1.*x2 + x1 + x2.^4 + x2.^3 + x2.^2 + x2);
% % u = K*(x1.^4 + (x1.^3).*x2 + x1.^3 + (x1.^2).*(x2.^2) + (x1.^2).*x2 + x1.^2 + x1.*(x2.^3) + x1.*(x2.^2) + x1.*x2 + x2.^4 + x2.^3 + x2.^2);
% dx1= - x1 + x2 + x1.^2 + x1.*x2 - x1.^3 + (x1.^2).*x2 - x1.*(x2.^2) + 2*x2.^3;
% dx2= - x1 + 1.5*x2 - x1.^2 - 0.5*x1.*x2 - x1.^3 - (x1.^2).*x2 + 0.5*x1.*(x2.^2) - 2*x2.^3 + u;

figure
streamslice(x1, x2, dx1, dx2, 1.5);
hold on
plot(0, 0, 'o')
hold on
plot(0, 0, 'x')
axis tight equal

return