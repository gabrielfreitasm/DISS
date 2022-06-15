clear
close all
clc

warning off

A = [-3, 1;
    -2, 0];

B = [0;
    1];

C = eye(size(A));

D = zeros(size(C,1),size(B,2));

% initial condition
x0 = 0;

P = sdpvar(size(A,1), size(A,2), 'symmetric', 'real');
Q = sdpvar(size(C,1), size(C,1), 'symmetric', 'real');
S = sdpvar(size(C,1), size(B,2), 'full', 'real');
R = sdpvar(size(B,2), size(B,2), 'symmetric', 'real');

LMI1 = P - 1e-5*eye(size(P)) >= 0;

LMI2 = [A'*P + P*A - C'*Q*C, P*B - C'*(Q*D + S);
    B'*P - (Q*D + S)'*C, - R - S'*D - D'*S - D'*Q*D] + 1e-5*eye(size(A,1) + size(R,1)) <= 0;

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
Q = value(Q);
S = value(S);
R = value(R);

return