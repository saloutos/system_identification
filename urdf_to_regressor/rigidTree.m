addpath('C:/Users/stare/casadi-windows-matlabR2016a-v3.5.5');
addpath(genpath('../ROAM-lab-inertia_id/'))
addpath(genpath('urdf2casadi'))
setPath;

import urdf2casadi.Utils.modelExtractionFunctions.*
% creating model structure for spatial_v2 (smds)
[smds, model] = extractSystemModel('SingleArm_URDF.urdf');

% creating rigidBodyTree
%robot = importrobot('SingleArm_URDF.urdf');

% Trying envs hack
% pyExec = 'C:\Users\stare\mambaforge\envs\python.exe';
% pyRoot = fileparts(pyExec);
% p = getenv('PATH');
% p = strsplit(p, ';');
% addToPath = {
%     pyRoot
%     fullfile(pyRoot, 'Library', 'mingw-w64', 'bin')
%     fullfile(pyRoot, 'Library', 'usr', 'bin')
%     fullfile(pyRoot, 'Library', 'bin')
%     fullfile(pyRoot, 'Scripts')
%     fullfile(pyRoot, 'bin')
%     };
% p = [addToPath(:); p(:)];
% p = unique(p, 'stable');
% p = strjoin(p, ';');
% setenv('PATH', p);