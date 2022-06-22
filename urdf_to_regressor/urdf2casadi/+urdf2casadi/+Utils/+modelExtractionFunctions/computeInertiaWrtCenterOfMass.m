function [I_oL, m, com ,rpy] = computeInertiaWrtCenterOfMass(model, jointIndex)
%Compute the inertia of body i wrt a frame centered in the center of mass.
% Seethe inertial tag for a link in a urdf:
% http://wiki.ros.org/urdf/XML/link
% Link 1 is the base_link, which is not considered in the array of
% spatial inertias(see http://royfeatherstone.org/spatial/v2/sysmodel.html)

% If a link does not have the inertial tag we'll set default values to zero
if ~isfield(model.robot.link{jointIndex+1}, 'inertial')
    warning('Setting inertial elements(center of mass,mass,Inertia Matix) to zero for link %d', jointIndex);
    com = [0 0 0];
    m = 0;
    I = zeros(3);
    rpy = [0 0 0];
else
    % Center of mass
    com = str2num(model.robot.link{jointIndex+1}.inertial.origin.Attributes.xyz);

    % Mass
    m = str2num(model.robot.link{jointIndex+1}.inertial.mass.Attributes.value);

    % Inertia
    ixx = str2num(model.robot.link{jointIndex+1}.inertial.inertia.Attributes.ixx);
    
    ixy = str2num(model.robot.link{jointIndex+1}.inertial.inertia.Attributes.ixy);
    
    ixz = str2num(model.robot.link{jointIndex+1}.inertial.inertia.Attributes.ixz);
    
    iyy  = str2num(model.robot.link{jointIndex+1}.inertial.inertia.Attributes.iyy);
    
    iyz = str2num(model.robot.link{jointIndex+1}.inertial.inertia.Attributes.iyz);
    
    izz = str2num(model.robot.link{jointIndex+1}.inertial.inertia.Attributes.izz);
    
    I =[ixx, ixy, ixz;
        ixy, iyy, iyz;
        ixz, iyz, izz]; 

    % Sanity check on Inertia: it must be at least positive semidefinite
    principalI = eig(I);
    if find(principalI<0)
        warning('Inertial matrix body %d is not positive semi-definite', jointIndex);
    end
    % Get the orientation of the frame (G) wrt the Inertia is computed in the urdf
    % and the local body frame
    if isfield(model.robot.link{1,jointIndex+1}.inertial.origin.Attributes,'rpy')
        rpy = str2num(model.robot.link{1,jointIndex+1}.inertial.origin.Attributes.rpy);
    else
        warning('Link %d missing orientation of the frame wrt the inertia is computed in the urdf.\n Setting it to [0 0 0]', jointIndex);
        rpy = [0 0 0];
    end
end
% Compute the Inertia matrix wrt a frame centered in the center of mass(g) but
% alligned with body local frame (L) orientation 
I_gL = urdf2casadi.Utils.modelExtractionFunctions.changeReferenceFrameInertiaMatrix(I, [],rpy,[], 'rpy');
% Compute Inertia matrix wrt frame centered in the local frame origin (o) and
% with orintation of the local body frame (L)
I_oL = urdf2casadi.Utils.Spatial.mcI(m, com, I_gL);   
end