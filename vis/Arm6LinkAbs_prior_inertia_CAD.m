function [pi_prior, P, Q_bound_ellipse] = Arm6LinkAbs_prior_inertia_CAD()

    % TODO: how to deal with bodies not having rotors? just make rotor small?

    [model, graphics]     = Arm6LinkAbsModel();
    N_RB = model.N_RB; % number of rigid bodies
    NJ = model.NB;
    NL = model.NL;
    
    pi_prior = zeros(10, N_RB);
    Q_bound_ellipse = cell(N_RB,1);
    P =cell(N_RB,1);    % Pseudo inertia matrix for body i
    Qc=cell(N_RB,1);    % 3x3 Matrix which provides a bounding ellipse for the full rigid body
    c =cell(N_RB,1);    % Center of Bounding Ellipsoid
    
    % Prior Params
    for i = 1:NL % for number of joints (should be N_RB/2 if all joints have rotors)
       % pi_prior orders all bodies then all rotors
       pi_prior(:,i)   = inertiaMatToVec(model.Il{i});
       pi_prior(:,i+NL) = inertiaMatToVec(model.Ir{i});
    end
    
    % Loop through one body at a time (including the rotors)
    for i = 1:N_RB
        P{i} = inertiaVecToPinertia(pi_prior(:,i));
        
        assert( min(eig(P{i})) > 0 )
        
        if i <= NL
            boundSemiAxisLengths   = graphics{i}.boundAxes;
            boundCenter            = graphics{i}.boundCenter;
        else
            boundSemiAxisLengths   = graphics{i-NL}.boundAxesRot;
            boundCenter            = graphics{i-NL}.boundCenterRot;
        end

        Q_bound_ellipse{i} = SemiAxesToQ( diag(boundSemiAxisLengths) ,boundCenter);

        % If the parameters a{i} are realizable with density on the bounding
        % ellipse, this trace must be positive
        tQP = trace(Q_bound_ellipse{i}*P{i});

        assert( tQP > 0 )
    end
end
