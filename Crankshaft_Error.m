function [sigma_v_mod] = Crankshaft_Error(r,d,l,theta,w,sigma_dist)
%CRANKSHAFT_ERROR Computes uncertainty in calculated model velocity
%
%   Inputs:     r                 -   distance between origin (rotation axis) and attachment point A [m]
%               d                 -   horizontal distance between vertical shaft and disk center [m]
%               l                 -   connecting bar length [m]
%               theta             -   disk angular position data [deg]
%               w                 -   disk angular velocity data [deg/s]
%               sigma_dist        -   uncertainty in the distance measurements [m]
%                                     r, d, and l [m] (all three have the same uncertainty)
%
%   Outputs:    sigma_v_mod       -   uncertainty in collar vertical speed [cm/s]
%

% Convert the angle theta from [deg] to [rad]
theta = deg2rad(theta);

% Convert angular velocity omega from [deg/s] to [rad/s]
w = deg2rad(w);

% BETA
% Calculate angle beta [rad]
beta = asin( (r - d .* sin(theta)) ./ l );

% beta is a function of r, d, and l, which all have associated uncertainty
partial_beta_r = ((1 - (((r-d.*sin(theta))./l).^2)).^(-0.5)) .* (-sin(theta)./l);
partial_beta_d = ((1 - (((r-d.*sin(theta))./l).^2)).^(-0.5))./l;
partial_beta_l = -((1 - (((r-d.*sin(theta))./l).^2)).^(-0.5)) .* ((r-d.*sin(theta))./l.^2);

% Each element that will be squared in the square root for the general method
e_r = partial_beta_r .* sigma_dist;
e_d = partial_beta_d .* sigma_dist;
e_l = partial_beta_l .* sigma_dist;

error_mat = [e_r, e_d, e_l]; % each row corresponds to one data point

sigma_beta = vecnorm(error_mat,2,2); % nx1 column vector containing the 2-norm (e.g. distance formula) of each row of the input matrix


% V MOD
% v_mod is a function of beta and r, which have their own associated uncertainty
partial_v_mod_r = w .* (cos(theta) .* tan(beta) - sin(theta));
partial_v_mod_beta = w .* d .* (cos(theta) .* (sec(beta).^2));

% Each element that will be squared in the square root for the general method
e_r = partial_v_mod_r .* sigma_dist;
e_beta = partial_v_mod_beta .* sigma_beta;

error_mat = [e_r, e_beta]; % each row corresponds to one data point

sigma_v_mod = vecnorm(error_mat,2,2); % nx1 column vector containing the 2-norm (e.g. distance formula) of each row of the input matrix
sigma_v_mod = 100 .* sigma_v_mod; % convert from [m/s] to [cm/s]
end

