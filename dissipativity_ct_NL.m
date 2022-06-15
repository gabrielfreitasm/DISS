clear
close all
clc

warning off

% system
% Example: McCourt and Antsaklis (continuous)

syms x1 x2 x3 u real

vars = [x1; x2; x3; u];
x = [x1; x2; x3];

f = [-(x1^3 + x1*x3^2)*(1 + x3^2);
    -(x1^2*x2 + x2)*(1 + x3^2);
    -(x3 + x1^2*x3)*(1 + x3^2) - 3*x3];

g = [0;
    0;
    1];

h = x3;

y = h;

prog = sosprogram(vars);

[prog, P] = sospolymatrixvar(prog, monomials(x, 0), [length(x), length(x)], 'symmetric');
[prog, Q] = sospolymatrixvar(prog, monomials(vars, 0), [length(h), length(h)], 'symmetric');
[prog, S] = sospolymatrixvar(prog, monomials(vars, 0), [length(h), length(u)]);
[prog, R] = sospolymatrixvar(prog, monomials(vars, 0), [length(u), length(u)], 'symmetric');

V = x'*P*x + 0.867*x1^4 + 0.815*x2^4 + 0*x3^4;

gradV = [diff(V, x1); diff(V, x2); diff(V, x3)];

diss = gradV'*(f + g*u) - y'*Q*y - 2*y'*S*u - u'*R*u;

prog = sosmatrixineq(prog, P - 1e-5*eye(size(P)), 'quadraticMineq');
prog = sosineq(prog, -diss);

sol = sossolve(prog);

V = sosgetsol(sol, V);
P = double(sosgetsol(sol, P));
Q = double(sosgetsol(sol, Q));
S = double(sosgetsol(sol, S));
R = double(sosgetsol(sol, R));

return