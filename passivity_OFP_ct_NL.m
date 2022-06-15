clear
close all
clc

warning off

% system
% Example: McCourt and Antsaklis (continuous)

syms x1 x2 x3 u real

vars = [x1; x2; u];
x = [x1; x2];

f = [-2*x1 + x2 - 0.5*x1^3;
    -0.5*x1 - x2 - x2^3];

g = [0;
    1];

h = x2 + 2*u;

y = h;

prog = sosprogram(vars);

[prog, V] = sospolyvar(prog, monomials(x, [2 4]));
[prog, rho] = sospolyvar(prog, monomials(x, 0));

gradV = [diff(V, x1); diff(V, x2)];

pass = gradV'*(f + g*u) + rho*y'*y - u'*y;

prog = sosineq(prog, V);
prog = sosineq(prog, -pass);

prog = sossetobj(prog, -rho);

sol = sossolve(prog);

V = sosgetsol(sol, V);
rho = double(sosgetsol(sol, rho));

return