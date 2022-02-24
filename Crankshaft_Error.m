function [sigma_Beta,sigma_VB] = Crankshaft_Error(d,r,L,theta,omega,Beta,sigma_d,sigma_r,sigma_L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
partial_Beta_d = ((1 - (((d-r*sin(theta))/L)^2))^(-1/2)) /L;
partial_Beta_r = ((1 - (((d-r*sin(theta))/L)^2))^(-1/2)) * (-sin(theta)/L);
partial_Beta_L = ((1 - (((d-r*sin(theta))/L)^2))^(-1/2)) * ((d-r*sin(theta))/L^2) * -1;
partial_VB_r = omega * (cos(theta)*tan(Beta) - sin(theta));
partial_VB_Beta = omega*r*(cos(theta)*(sec(Beta)^2));
error_Beta_d = partial_Beta_d * sigma_d;
error_Beta_r = partial_Beta_r * sigma_r;
error_Beta_L = partial_Beta_L * sigma_L;

 sigma_Beta = norm([error_Beta_d error_Beta_r error_Beta_L]);
error_VB_r = partial_VB_r * sigma_r;
error_VB_Beta = partial_VB_Beta * sigma_Beta;
 sigma_VB = norm([error_VB_r error_VB_Beta]);


end

