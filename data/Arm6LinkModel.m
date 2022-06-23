function [model, graphics] = Arm6LinkModel()

% This function creates both the rigid body model struct and a graphics
% cell array that specifies the full teleop arm model

% The 0-configuration for the robot is completely vertical with all joint
% frames / body frames aligned
% The inertial coordinates have +z up (i.e. gravity = [0 0 -9.81 m/s^2])
% The body coordinates have +x forward, +y left, +z up

model.gc.point = zeros(3,0); % what do these mean?
model.gc.body = zeros(1,0);

% There are 6 links/bodies...12 total bodies with the WP and WR rotors, 10 without
Nb = 0; % start body index at 0

%% Nominal arm parameters

% motor sizes
U12_dia = 0.120;
U10_dia = 0.100;
U12_width = 0.04;
U10_width = 0.035;

% masses (approx from URDF)
link_masses = [2.93, 1.74, 4.99, 0.82, 0.17, 0.65]; % added 600g for gripper to last mass
rotor_masses = [0.1, 0.06]; % update these from CAD?

% link lengths/offsets (from URDF) (first is at origin, then from prev. joint origin)
joint_offsets = [0, 0, 0; %0, 0, 0.051;
                 0, 0, 0.106;
                 0, 0, 0.071;
                 0, -0.0095, 0.3855;
                 0, 0.004, 0.035;
                 0, 0, 0];

% rotor offsets (from CAD) (from corresponding joint origin for joints 1-4, from upper arm origin for 5 and 6)
rotor_offsets = [0, 0, -0.02125;
                0, 0.09625, 0;
                0, 0, -0.02125;
                0, 0.07088, -0.307;
                0, 0, 0; % not using these rotors for now
                0, 0, 0];

% link centers (approx from CAD) used for bounding boxes, inertia estimates
% TODO: these might not be necessary?
link_centers = [0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0;
                0, 0, 0];

% bounding boxes (approx from CAD)
base_box = [0.125, 0.2, 0.175]; % all box lengths are [x,y,z] % TODO: add center locations here?
shoulder_box = [0.130, 0.138, 0.140];
upper_arm_box = [0.130, 0.140, 0.430];
lower_arm_box = [0.070, 0.050, 0.420];
wrist_box = [0.040, 0.060, 0.060];
gripper_box = [0.050, 0.110, 0.200]; % represents wrist roll output

%% Shoulder 1 (Shoulder yaw output link)

Nb = Nb+1;
model.parent(Nb) = Nb-1;
model.jtype{Nb}  = 'Rz';
model.jtype_rotor{Nb}  = 'Rz';
model.Xtree{Nb}  = plux(eye(3), joint_offsets(1,:)');
model.Xrotor{Nb} = plux(eye(3), rotor_offsets(1,:)');

model.I{Nb}      = mcI( masses(1), joint_offsets(2,:)/2, boxInertia(masses(1),base_box) );
model.I_rotor{Nb}   = mcI( rotor_masses(1), [0 0 0], boxInertia(rotor_masses(1), [U12_dia, U12_dia, U12_width]) );

model.gr{Nb} = 6;

graphics{Nb}.boundCenter = joint_offsets(2,:)'/2;
graphics{Nb}.boundAxes   = base_box/2*1.8;
graphics{Nb}.boundCenterRot = [0 0 0]';
graphics{Nb}.boundAxesRot   = ([U12_dia, U12_dia, U12_width]/2)*1.8;

%% Shoulder 2 (Shoulder pitch output link)

%% Upper arm (Shoulder roll output link)

%% Lower arm (Elbow output link)

%% Wrist (Wrist pitch output link)

%% Gripper (Wrist roll output link)

% adding rotor to this one will be extra tricky...need to rotate output
% frame to get positive axis?



%% Wrap up
mT = 0; % what is this for?
for i = 1:length(model.I)
    mT = mT+model.I{i}(6,6);
end
%mT

model.NB   = Nb;

model.has_rotor=[1 1 1 1 0 0]'; % TODO: ignoring wrist rotors for now!
model = postProcessModel(model);

end

function I = boxInertia(mass, x) % inertia for modeling as box of mass m and lengths x with uniform density
    I = (norm(x)^2*eye(3) - diag(x.^2))*mass/12;
end