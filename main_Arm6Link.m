% Need to run this file with its enclosing folder as the working directory
fileInfo = dir(matlab.desktop.editor.getActiveFilename);
cd(fileInfo.folder);

clear
close all
addpath(genpath([pwd,'/data']))
addpath(genpath([pwd,'/utility']))
addpath(genpath([pwd,'/vis']))

%% Load Data

% TODO: load correct data here, make sure formats are correct too
load('armlog2_processed_data.mat');

model = Arm6LinkModel();

%% Load Object files
disp("loading CAD mesh files (this may take a while)")
obj = cell(7,1);

obj{1} = readObj('Link1_ShoulderYaw.obj');
obj{2} = readObj('Link2_ShoulderPitch.obj');
obj{3} = readObj('Link3_ShoulderRoll.obj');
obj{4} = readObj('Link4_Elbow.obj');
obj{5} = readObj('Link5_WristPitch.obj');
obj{6} = readObj('Link6_WristRoll_WithGripper.obj');
obj{7} = readObj('Link0_BaseLink.obj');

disp("loaded CAD mesh files")

%% Data Setup
N = length(q);
n_links = 6; % # of conventional links
n_rotors= 6; % # of rigid-bodies modeled for rotors
n_bodies= n_links+n_rotors;
n_dofs  = length(qd{1});

% Time vector is saved in dataset!

% Need to generate the regressors since they are not in our dataset
regenerate_regressors = 1;
if regenerate_regressors
   Y = {};
   for i =1:N
       if mod(i,100) == 0
           fprintf('%d / %d Regressors Computed\n',i,N);
       end
       [Y_i, Yrot_i] = RegressorClassical( model, q{i}, qd{i}, qdd{i});
       Y{i,1} = [Y_i Yrot_i];
    end
end

% Setup spring/counterbalance regressor (for 6-link arm model)
Ks = repmat({zeros(n_dofs,1)}, N,1);
for i = 1:N
    Ks{i} = [0; cb_torque_coef(q{i}); 0; 0; 0; 0];
end

% Setup Friction Regressor
B = repmat({zeros(n_dofs,n_dofs)}, N,1 );
Bc = repmat({zeros(n_dofs,n_dofs)}, N,1 );
for i = 1:N
   B{i}  = diag( qd{i} );
   Bc{i} = diag( sign( qd{i} ) );
end

% Convert cell arrays to matrices
tau_stack = cell2mat(tau);
Y_stack = cell2mat(Y);
tau_mat = cell2mat(tau')';
B_stack = cell2mat(B);
Bc_stack = cell2mat(Bc);
Ks_stack = cell2mat(Ks);

% TODO: include spring regressor here
Y_total = [Y_stack B_stack Bc_stack Ks_stack];

% Optional: Select subset for training data
N_train   = N; % Full data set
Y_train   =   Y_stack(1:n_dofs*N_train,:);
B_train   =   B_stack(1:n_dofs*N_train,:);
Bc_train  =  Bc_stack(1:n_dofs*N_train,:);
tau_train = tau_stack(1:n_dofs*N_train,:);
Ks_train = Ks_stack(1:n_dofs*N_train,:);

%% Prior
% The parameter order for each link is: 
% m, h_x, h_y, h_z, I_xx, I_yy, I_zz, I_yz, I_xz, I_xy

[pi_prior, ... % Prior parameters
    J_prior, ... % Prior Pseudo-inertias
      Q_bound] ... % Bounding Ellipsoids S={ x | [x ; 1]'*Q*[x ; 1] >= 0 }
      = Arm6Link_prior_inertia_CAD();

figure(1); clf;
color = rand(12,3);
title('Bounding Ellipsoids');
Arm6Link_bound_visualize(obj, Q_bound, color)
grid on;
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');

  
%% Conventional System ID with Entropic Regularization

fprintf('=========== Convex Approach (SDP with LMIs) =========\n');

cvx_setup;

% Clear CVX
cvx_begin
cvx_end

% Set solver
cvx_solver mosek % <- Can be changed if you don't have it.

weight_regularization = 1e-1;
use_const_pullback_approx = 0;
force_ellipsoid_density_realizability = 1;

cvx_begin
    variable params(10,n_bodies)  % inertial params of links / rotors (units vary)
    variable b(n_dofs)            % viscous friction coefficient (Nm / (rad/s) 
    variable bc(n_dofs)           % coulomb friction coefficient (Nm) 
    % TODO: add variable for spring constant...does this initialize a scalar?
    variable ks                   % efffective spring contstance (N/mm)
    
    variable J(4,4,n_bodies) semidefinite % pseudo-inertia for each body
    
    expression bregman(n_bodies) % Entropic / bregman divergence for regularization
    expression e_train(3*N_train,1); % Training error (Nm)
    
    % Training error
    e_train = Y_train*params(:) + B_train*b + Bc_train*bc + Ks_train*ks - tau_train;
    % TODO: add variable for spring constant

    % Entropic (i.e., bregman) divergence between psuedo-inertias
    for i=1:n_bodies
        if use_const_pullback_approx == 0
            % d( J || J_prior )
            bregman(i) = -log_det( J(:,:,i) ) + log(det(J_prior{i})) ...
                            + trace(J_prior{i} \ J(:,:,i) ) - 4;
        else
            % constant pullback approximation
            M{i} = pullbackMetric(pi_prior(:,i));
            bregman(i) = 1/2*(params(:,i) - prior_params(:,i))'*M{i}*(params(:,i) - prior_params(:,i));
        end
    end
    
    % Objective = Squared 2-norm of residual + Regularization
    minimize( 1/2*sum_square_abs(e_train)/length(e_train) ...
                    + weight_regularization*sum(bregman) )
    
    % Only constraints are to set the pseudo-inertias based on params
    subject to:
        for i=1:n_bodies
            J(:,:,i) == inertiaVecToPinertia( params(:,i) );
            if force_ellipsoid_density_realizability
                trace( J(:,:,i)*Q_bound{i}  ) >= 0;
            end
        end
        % TODO: add constraint for positive spring constant, positive friction coefficients?
        ks >= 0;
        

cvx_end

%%
tau_predict_entropic = reshape(Y_stack*params(:) + B_stack*b + Bc_stack*bc + Ks_stack*ks, n_dofs, N)';
plotTorquePredictions(2,'Convex, Entropic Regularized',t,tau_mat, tau_predict_entropic);
% TODO: add variable for spring constant

% Visualize inertia and CAD
color = rand(12,3);
figure(3); clf;

% Represents the inertial params by solid ellipsoid of uniform density
% Note: This solid need not be within the bounding ellipsoid since
%       it represents but one density distribution that realizes the 
%       parameters
Arm6Link_ellipsoid_visualize(obj, params, color);

grid on;
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
title('Inertia Ellipsoids (Convex, Entropic Regularized)');


%% Kinematics Plotting
q_mat = cell2mat(q')';
figure(6); clf;
subplot(321);
plot(t,q_mat(:,1))
ylabel('Sh. yaw angle (rad)');
subplot(323);
plot(t,q_mat(:,2))
ylabel('Sh. pitch angle (rad)');
subplot(325);
plot(t,q_mat(:,3))
ylabel('Sh. roll angle (rad)');
xlabel('Time (s)');

subplot(322);
plot(t,q_mat(:,4))
ylabel('Elbow angle (rad)');
subplot(324);
plot(t,q_mat(:,5))
ylabel('Wrist pitch angle (rad)');
subplot(326);
plot(t,q_mat(:,6))
ylabel('Wrist roll angle (rad)');
xlabel('Time (s)');


%% Helpers
function plotTorquePredictions(figNum, name, t, tau_actual , tau_predict)
    tau_rms = rms(tau_actual - tau_predict);
    figure(figNum)
    clf
    subplot(321)
    plot(t,tau_actual(:,1)); hold on;
    plot(t,tau_predict(:,1),'r','LineWidth',1.5 )
    ylabel('Sh. Yaw Torque (Nm)');
    rms_tag = sprintf('(%s) [RMS=%.2f, %.2f, %.2f (Nm)]',name, tau_rms(1), tau_rms(2), tau_rms(3));
    title(['Torque Predictions ' rms_tag]);
    legend('Measured','Predicted')

    subplot(323)
    plot(t,tau_actual(:,2)); hold on;
    plot(t,tau_predict(:,2),'r','LineWidth',1.5)
    ylabel('Sh. Pitch Torque (Nm)');

    subplot(325)
    plot(t,tau_actual(:,3)); hold on;
    plot(t,tau_predict(:,3),'r','LineWidth',1.5)
    ylabel('Sh. Roll Torque (Nm)');
    xlabel('Time (s)');

    subplot(322)
    plot(t,tau_actual(:,4)); hold on;
    plot(t,tau_predict(:,4),'r','LineWidth',1.5)
    ylabel('Elbow Torque (Nm)');
    xlabel('Time (s)');

    subplot(324)
    plot(t,tau_actual(:,5)); hold on;
    plot(t,tau_predict(:,5),'r','LineWidth',1.5)
    ylabel('Wrist Pitch Torque (Nm)');
    xlabel('Time (s)');

    subplot(326)
    plot(t,tau_actual(:,6)); hold on;
    plot(t,tau_predict(:,6),'r','LineWidth',1.5)
    ylabel('Wrist Roll Torque (Nm)');
    xlabel('Time (s)');
end

function x = cb_torque_coef(q)
    % takes vector of joint angles, returns the shoulder pitch counterbalance torque divided by the spring constant
    joint_angle = q(2); 
    a = 11; %arm attachment distance
    l = 80; %spring free length
    b = -52; %pivot attachment distance
    x_s1 = -a*sin(joint_angle);
    y_s1 = a*cos(joint_angle);
    s1 = sqrt(x_s1^2 + (y_s1-b)^2);
    spring_angle = acos(a^2 + s1^2 - b^2)/(2*a*s1);
    x = (l-s1) * a * 0.001 * sin(spring_angle);
    if (joint_angle>0)
        x = -x;
    end
    % Note: tau_spring = ks*x
end



