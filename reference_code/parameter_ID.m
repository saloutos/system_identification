function sol = parameter_ID(block,params_prior,Y,tau_vec,c_true,opt_type)
%parameter_ID   Formulates and solves an optimization problem for
%identifying inertial parameters of for the block simulation

%   no friction right now
%   inputs: 
%           true inertia parameters, block
%           initial guess at inertia parameters, params_prior
%           regressor matrix, Y
%           force vector, tau_vec
%           opt_type, 'LMI','ED',or'CP' to determine optimization
%   outputs:
%           optimization result, sol
    
    verbosity = 0;

    % LMI optimization

    if strcmp(opt_type,'LMI')

        inertia_params = sdpvar(10,1); % m, m*c (3x1), I(6x1)
        eps_J = sdpvar(1);
        % Note: I is [Ixx, Ixy, Ixz, Iyy, Iyz, Izz], taken at point of applied forces

        params = [inertia_params];
        gamma = 0.1;

        % res = [Yreg,Vreg]*params - tau_vec; % with friction
        res = Y*params - tau_vec; % without friction
        res_params = params - params_prior;
        LMI_objective = res'*res + eps_J*eps_J;
        LMI_objective = LMI_objective + gamma*(res_params'*res_params); % add regularization

        % physically consistent constraints
        Qs = diag([100,100,100]);
        Qs_4 = diag([100,100,100,1]);
        C = [params(1), params(2:4)'; params(2:4), Qs];
        Ibar = Ivect_to_mat(params(5:10));
        J = [0.5*trace(Ibar)*eye(3) - Ibar, params(2:4); params(2:4)', params(1)];

        % LMI_constraints = [0<=eps_J<=1e-3, C>=0, J>=eps_J, trace(J*Qs_4)>=0];

        LMI_constraints = [0<=eps_J<=1e-3, J>=eps_J];

        ops = sdpsettings('verbose',verbosity);

        optimize(LMI_constraints,LMI_objective,ops);
        sol = value(params);
        mhat = sol(1);
        c_hat = sol(2:4)/sol(1);
        Ip = Ivect_to_mat(sol(5:10));
        Ic_hat = Ip - mhat*(c_hat'*c_hat*eye(3) - c_hat*c_hat');

        % display some info

%         fprintf('____LMI Optimization Results____\n\n');
% 
%         disp('Slack:')
%         disp(value(eps_J))
% 
%         disp('Masses:')
%         disp(mhat); disp(block.M); 
% 
%         disp('CoMs:')
%         disp([c_hat, c_true']);
% 
%         disp('Inertia:')
% %         disp(Ip); 
%         disp(Ic_hat); disp(block.Ic); 

    %     disp('Friction:')
    %     disp(mu_vec)
    
    
    % Entropic Divergence optimization

    % what is difference between regularized or point to set formulations?
    % using regularized formulation

    elseif strcmp(opt_type,'EntDiv')

        % start with no friction
        params = sdpvar(10,1); % m, m*c (3x1), I(6x1)
        gamma = 0.1;
        % Note: I is [Ixx, Ixy, Ixz, Iyy, Iyz, Izz], taken at point of applied forces

        sig = 1; %*eye(length(tau_vec)); % observation variance...should be Nsamp x Nsamp

        % only one body
        P = [0.5*trace(Ivect_to_mat(params(5:10)))*eye(3) - Ivect_to_mat(params(5:10)), params(2:4); params(2:4)', params(1)];
        P_prior = [0.5*trace(Ivect_to_mat(params_prior(5:10)))*eye(3) - Ivect_to_mat(params_prior(5:10)), params_prior(2:4); params_prior(2:4)', params_prior(1)];

        A = sig^(-0.5) * Y;
        b = sig^(-0.5) * tau_vec;

        res = A*params - b;

        J_ls = res'*res;
        % regularizer is entropic divergence
        J_reg = -geomean(P) + trace(P_prior\P);
        ED_objective = J_ls + gamma*J_reg;

        % can add constraints for E-density realizability
        Qs_4 = diag([100,100,100,1]);
%         ED_constraints = [];
        ED_constraints = [trace(Qs_4*P)>= 0];

        ops = sdpsettings('verbose',verbosity);
        optimize(ED_constraints, ED_objective, ops);

        sol = value(params);
        mhat = sol(1);
        c_hat = sol(2:4)/sol(1);
        Ip = Ivect_to_mat(sol(5:10));
        Ic_hat = Ip - mhat*(c_hat'*c_hat*eye(3) - c_hat*c_hat');

        % display some info

%         fprintf('____Entropic Divergence Optimization Results____\n\n');
% 
%         disp('Masses:')
%         disp(mhat); disp(block.M); 
% 
%         disp('CoMs:')
%         disp([c_hat, c_true']);
% 
%         disp('Inertia:')
% %         disp(Ip); 
%         disp(Ic_hat); disp(block.Ic); 

    

    % Constant Pullback optimization

    % using regularized formulation
    % add physical constistency constraints?
    
    elseif strcmp(opt_type,'ConPb')

        params = sdpvar(10,1); % m, m*c (3x1), I(6x1)
        gamma = 0.1;
        % Note: I is [Ixx, Ixy, Ixz, Iyy, Iyz, Izz], taken at point of applied forces

        sig = 1; %*eye(length(tau_vec)); % observation variance...should be Nsamp x Nsamp

        % only one body
        P = [0.5*trace(Ivect_to_mat(params(5:10)))*eye(3) - Ivect_to_mat(params(5:10)), params(2:4); params(2:4)', params(1)];
        P_prior = [0.5*trace(Ivect_to_mat(params_prior(5:10)))*eye(3) - Ivect_to_mat(params_prior(5:10)), params_prior(2:4); params_prior(2:4)', params_prior(1)];

        A = sig^(-0.5) * Y;
        b = sig^(-0.5) * tau_vec;

        const_pbmet_sqrt = pullback_metric(P_prior)^(0.5);

        res = A*params-b;
        J_ls = res'*res;
        % constant pullback regularizer
        res_pb = (params - params_prior)' * const_pbmet_sqrt * (params - params_prior);
        J_reg = res_pb'*res_pb;
        CP_objective = J_ls + gamma*J_reg;

        % can add constraints for E-density realizability
        Qs_4 = diag([100,100,100,1]);
        CP_constraints = [trace(Qs_4*P)>= 0];
        % TODO: constraint on pseudo-inertias as well?

        ops = sdpsettings('verbose',verbosity);
        optimize(CP_constraints, CP_objective,ops);
        
        sol = value(params);
        mhat = sol(1);
        c_hat = sol(2:4)/sol(1);
        Ip = Ivect_to_mat(sol(5:10));
        Ic_hat = Ip - mhat*(c_hat'*c_hat*eye(3) - c_hat*c_hat');
        
        % display some info

%         fprintf('____Constant Pullback Optimization Results____\n\n');
% 
%         disp('Masses:')
%         disp(mhat); disp(block.M); 
% 
%         disp('CoMs:')
%         disp([c_hat, c_true']);
% 
%         disp('Inertia:')
% %         disp(Ip); 
%         disp(Ic_hat); disp(block.Ic); 

    else
       fprintf('Invalid type of optimization formulation. Try LMI, EntDiv, or ConPb.\n'); 
    end

end

function g = pullback_metric(P)
    % takes pseudo-inertia P and returns pullback metric matrix
    g = zeros(10,10);
    P_inv = inv(P);
    for ii=1:10
        for jj=1:10
            e_i = zeros(10,1); e_i(ii) = 1;
            e_j = zeros(10,1); e_j(jj) = 1;
            V_i = [0.5*trace(Ivect_to_mat(e_i(5:10)))*eye(3) - Ivect_to_mat(e_i(5:10)), e_i(2:4); e_i(2:4)', e_i(1)];
            V_j = [0.5*trace(Ivect_to_mat(e_j(5:10)))*eye(3) - Ivect_to_mat(e_j(5:10)), e_j(2:4); e_j(2:4)', e_j(1)];
            g(ii,jj) = trace(P_inv * V_i * P_inv * V_j);     
        end
    end    
end

function Imat = Ivect_to_mat(Ivect)
    Imat = [Ivect(1), Ivect(2), Ivect(3); Ivect(2), Ivect(4), Ivect(5); Ivect(3), Ivect(5), Ivect(6)];
end

