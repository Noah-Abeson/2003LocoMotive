%% Header
% Last Modified: 2/22/22

clear; close all;
%% Main Body

% Manipulate disk data

num_files = 6;
files = ["Test1_10pt5V","Test1_9pt5V","Test1_8pt5V","Test1_7pt5V","Test1_6pt5V","Test1_5pt5V"];

for ii = 1:num_files
    [theta_exp{ii},w_exp{ii},v_exp{ii},t_exp{ii}] = LCSDATA(files(ii));
end

%% Plotting
for ii = 1:num_files
    figure(); hold on;
    plot(theta_exp{ii},v_exp{ii});
    xlabel("Theta \theta [°]");
    ylabel("Velocity [cm/s]");
    title(files{ii},"Interpreter","none");
end
%% Function: LCSDATA
function [theta_exp,w_exp,v_exp,t_exp] = LCSDATA(filename)
%LCSDATA Reads in raw data file and outputs manipulated angular position and velocity of disk and vertical speed of collar.
%
%   Inputs:     filename          -   name of data file to import
%
%   Outputs:    theta_exp         -   experimental disk angular position data [deg]
%               w_exp             -   experimental disk angular velocity data [deg/s]
%               v_exp             -   experimental collar vertical speed data [cm/s]
%

% Read Data
data = readmatrix(filename);

% Extract Data
time = data(:,1);
angle = data(:,2);
wheel_speed = data(:,4);
slide_speed = data(:,5);

% Start theta between 0 and 360
start_angle = mod(angle(1),360);
angle_diff = angle(1) - start_angle;
theta_exp = angle - angle_diff;

% Find 6 revolutions
logic_vec = theta_exp < (6 * 360 + start_angle);

% Return values
theta_exp = theta_exp(logic_vec);
w_exp = wheel_speed(logic_vec);
v_exp = slide_speed(logic_vec) ./ 10;
t_exp = time(logic_vec);


end