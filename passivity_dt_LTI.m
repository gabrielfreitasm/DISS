clear
close all
clc

warning off

% system
% Example: McCourt and Antsaklis (discrete)
A = [0, 1;
    -2, -2];

B = [0;
    1];

C = [-1, 2];

D = 1.5;

dsys = ss(c2d(ss(A, B, C, D), 1, 'zoh'));
A = dsys.A;
B = dsys.B;
C = dsys.C;
D = dsys.D;

% initial condition
x0 = 0;

P = sdpvar(size(A,1), size(A,2), 'symmetric', 'real');

LMI1 = P - 1e-5*eye(size(P)) >= 0;

LMI2 = [A'*P*A - P, A'*P*B - C';
    B'*P*A - C, - D - D'] <= 0;

LMIs = [LMI1, LMI2];

options = sdpsettings('solver','sedumi','verbose',1);
sol = optimize(LMIs, [], options);

if sol.problem == 0
    [primal,~] = check(LMIs);
    if (min(primal) >= 0 && all(primal(1) > 0))
        disp('Sucessfully solved LMIs without problems'); 
    else
        disp('LMIs not solved');
    end
else
    [primal,~] = check(LMIs);
    if (min(primal) >= 0 && all(primal(1) > 0))
        disp(['Sucessfully solved LMIs, but solver acused ' yalmiperror(sol.problem)]);
    else
        disp(['LMIs not solved. Solver acused ' yalmiperror(sol.problem)]);
    end
end

P = value(P);

return