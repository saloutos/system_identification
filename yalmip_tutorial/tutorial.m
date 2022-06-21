% 4/1/19
% practice with yalmip (with mosek) tutorials

%% run yalmip test

yalmiptest

%% "Getting Started" tutorial

% Define variables
x = sdpvar(10,1);

% Define constraints 
Constraints = [sum(x) <= 10, x(1) == 0, 0.5 <= x(2) <= 1.5];
for i = 1 : 7
  Constraints = [Constraints, x(i) + x(i+1) <= x(i+2) + x(i+3)];
end

% Define an objective
Objective = x'*x+norm(x,1);

% Set some options for YALMIP and solver
% options = sdpsettings('verbose',1,'solver','quadprog','quadprog.maxiter',100);
options = sdpsettings('verbose',1,'solver','mosek'); %,'quadprog.maxiter',100);

% Solve the problem
sol = optimize(Constraints,Objective,options);

% Analyze error flags
if sol.problem == 0
 % Extract and display value
 solution = value(x)
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end

%% Semidefinite programming

A = [-1 2 0; -3 -4 1; 0 0 -2];
P = sdpvar(3,3);

F = [P >= 0, A'*P + P*A <=0];
F = [F, trace(P)==1];

% optimize(F);
% Pfeasible = value(P);
% check(F);

F = [F, P([2 3 6])>=0];
optimize(F,P(1,1));
value(P)

%% Determinant maximization example

n = 2;
P1 = randn(2,2); P1 = P1*P1'; % generate 2 random ellipsoids
P2 = randn(2,2); P2 = P2*P2'; 

t = sdpvar(2,1);
P = sdpvar(n,n);

F = [ 1>=t(1)>=0, 1>=t(2)>=0, (t(1)*P1)-P>=0, (t(2)*P2)-P>=0 ];
sol = optimize(F,-geomean(P));

x = sdpvar(2,1);
figure;
plot(x'*value(P)*x <= 1,[],'b'); hold on
plot(x'*value(P1)*x <= 1,[],'r');
plot(x'*value(P2)*x <= 1,[],'y');

%% Quadratic programming example

% linear regression

x = [1 2 3 4 5 6]';
t = (0:0.02:2*pi)'; 
A = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
e = (-4 + 8*rand(length(t),1)); % add noise
e(100:115) = 30; % add outliers

y = A*x + e; % noisy output data
%figure; plot(t,y)

% estimate vector x
xhat = sdpvar(6,1);
res = y - A*xhat;

% 1-norm (abs. val. of residuals)
bound = sdpvar(length(res),1);
Constraints = [-bound <= res <= bound];
optimize(Constraints,sum(bound)); % objective is bound of residuals
x_L1 = value(xhat);

% 2-norm (least-squares)
optimize([],res'*res); % QP with no constraints
x_L2 = value(xhat);

% inf-norm (minimize largest residual)
bound = sdpvar(1,1);
Constraints = [-bound <= res <= bound];
optimize(Constraints,bound);
x_Linf = value(xhat);

% plot all solutions
figure; plot(t, [y, A*x_L1, A*x_L2, A*x_Linf]);
legend('Noisy','L1','L2','Linf');




