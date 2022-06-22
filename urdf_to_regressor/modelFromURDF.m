function model = modelFromURDF(urdf)

% Takes urdf file, parses it into a rigidBodyTree, and then constructs a
% model

    addpath(genpath('../ROAM-lab-inertia_id/'))
    robot = importrobot(urdf);

    persistent memory;
    if length(memory) ~= 0
      model = memory;
      return
    end

    model.NB = robot.NumBodies;
    model.parent = [0 1];
    
    model.jtype = {'R', 'Ry'};	% 1st joint will be replaced by floatbase
    model.Xtree = {eye(6), eye(6)};
    
    model.gravity = [0 0 0];		% zero gravity is not the default,
                                            % so it must be stated explicitly
    
    % centroidal inertia of a unit-mass box with dimensions 0.2, 0.1 and 1:
    
    Ic = diag( [1.01 1.04 0.05] / 12 );
    
    model.I = { mcI(1,[0 -0.05 0],Ic),  mcI(1,[0 0.05 0],Ic) };
    
    % Draw the bodies
    
    model.appearance.body{1} = ...
      { 'box', [-0.1 -0.1 -0.5; 0.1 0 0.5], ...
        'colour', [0.6 0.3 0.8], ...
        'cyl', [0 -0.13 0; 0 0.13 0], 0.07 };
    
    model.appearance.body{2} = ...
      { 'box', [-0.1 0 -0.5; 0.1 0.1 0.5] };
    
    model.camera.direction = [0 -3 1];	% a better viewing angle
    model.camera.zoom = 0.9;
    
    model = floatbase(model);		% replace joint 1 with a chain of 6
                                            % joints emulating a floating base
    memory = model;

end

