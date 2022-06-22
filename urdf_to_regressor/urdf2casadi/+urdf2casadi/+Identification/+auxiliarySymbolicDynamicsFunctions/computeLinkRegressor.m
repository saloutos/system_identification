function A = computeLinkRegressor(linkSpatialVel, linkSpatialAcc)
% Compute the regressor associated with a single link. The regressor is the
% matrix that renders the dynamics Newton-Euler equation linear in the
% inertia parameters

% Import necessary functions
import urdf2casadi.Identification.auxiliarySymbolicDynamicsFunctions.computeRotationalMomentumRegressor
import urdf2casadi.Utils.Spatial.skew

% Link angular vel
linkAngularVel = linkSpatialVel(1:3);
% Link origin linear vel
linkLinearVel =  linkSpatialVel(4:6);
% Derivative of the components of the agular vel
linkAngularVelDot = linkSpatialAcc(1:3);
% Derivative of the velocity prat in the spatial acceleration
linkLinearVelDot = linkSpatialAcc(4:6);
% L(linkAngularVel)vec(Inertial_i) = Inertia_i*linkAngularVel
% Note that the Inertia must be expressed in the body i local frame
L_omega = computeRotationalMomentumRegressor(linkAngularVel);
L_omegaDot = computeRotationalMomentumRegressor(linkAngularVelDot);

% Compute the classical acceleration of the origin of body i frame(d_oi)
d_oi_DDot = linkLinearVelDot + skew(linkAngularVel)*linkLinearVel;

A = [zeros(3,1),-skew(d_oi_DDot),                    L_omegaDot+skew(linkAngularVel)*L_omega;
     d_oi_DDot, skew(linkAngularVelDot)+skew(linkAngularVel)*skew(linkAngularVel), zeros(3,6)];
end