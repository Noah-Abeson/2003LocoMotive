%% Header
% Last Modified: 2/22/22

clear; close all;
%% Data
% Experimental constants
r = 0.075; % [m]
d = 0.155; % [m]
l = 0.26; % [m]
sigma = 0.0005; % [m]

% Names of data files
num_files = 6;
files = ["Test1_10pt5V","Test1_9pt5V","Test1_8pt5V","Test1_7pt5V","Test1_6pt5V","Test1_5pt5V"];

% Read in and extract the experimental data: theta [deg], omega [deg/s],velocity [cm/s], time [s]
for ii = 1:num_files
    [theta_exp{ii},w_exp{ii},v_exp{ii},t_exp{ii}] = LCSDATA(files(ii));
end

%% Experimental Model
% Compute vertical velocity prediction: v_model [cm/s]
for ii = 1:num_files
    v_mod{ii} = LCSMODEL(r,d,l,theta_exp{ii},w_exp{ii});
end

%% Residuals
for ii = 1:num_files
    % Calculate the residuals (the difference between the experimental data
    % and the data predicted from the model
    res = v_exp{ii} - v_mod{ii};
    residual{ii} = res;
    
    % Calculate the mean and standard deviation of the residuals to
    % determine outliers
    avg = mean(res);
    average{ii} = avg; % store in cell array

    std_dev = std(res);
    standard_deviation{ii} = std_dev; % store in cell array

    % Remove the outliers
    non_outliers_logic = (res - avg) < (2 * std_dev);
    res_no_outlier = res(non_outliers_logic);
    residual_no_outlier{ii} = res_no_outlier; % store in cell array
    
    % Recalculate mean and standard deviation without outliers
    average_no_outlier{ii} = mean(res_no_outlier);
    standard_deviation_no_outlier{ii} = std(res_no_outlier); % store in cell array

    % Uncertainty
    sigma_v_mod{ii} = Crankshaft_Error(d,r,l,theta_exp{ii},w_exp{ii},sigma);
end

%% Plotting
for ii = 1:num_files
    
    % Vertical Velocity
    figure(); hold on; grid minor;
    plot(theta_exp{ii},v_exp{ii}); % Experimental
    plot(theta_exp{ii},v_mod{ii}); % Model
    xlabel("Angular Position \theta [°]");
    ylabel("Vertical Velocity [cm/s]");

    % Convert filename to voltage string for the title
    filename = files{ii};
    fprintf("%s\n",filename);
    voltage = filename(7:9);
    if voltage(2) ~= 'p'
        integer = str2num(voltage(1:2));
    else
        integer = str2num(voltage(1));
    end
    voltage_str = num2str(integer + 0.5);

    title(['Velocity vs Angle: ',voltage_str,'V']);
    legend("Experimental", "Model");

    % Residuals
    figure(); grid minor;
    plot(t_exp{ii},residual{ii})
    xlabel("Time [s]");
    ylabel("Differential Vertical Velocity [cm/s]");
    title(['Residual vs Time: ',voltage_str,'V']);

    % Error Bar
    figure(); grid minor;
    errorbar(t_exp{ii},residual{ii},sigma_v_mod{ii});
    xlabel("Time [s]");
    ylabel("Differential Vertical Velocity [cm/s]");
    title(['Residual vs Time with Uncertainty: ',voltage_str,'V']);

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
%% Function: LCSMODEL
function [v_mod] = LCSMODEL(r, d, l, theta, w)
%LCSMODEL Computes vertical speed of collar from locomotive crankshaft geometric constants and angular position and velocity of disk.
%
%   Inputs:     r                 -   distance between origin (rotation axis) and attachment point A [m]
%               d                 -   horizontal distance between vertical shaft and disk center [m]
%               l                 -   connecting bar length [m]
%               theta             -   disk angular position data [deg]
%               w                 -   disk angular velocity data [deg/s]
%
%   Outputs:    v_mod             -   collar vertical speed [cm/s]
%

% Calculate angle of arm, beta in degrees
beta = asind( (d - r .* sind(theta)) ./ l );
% Calculate velocity of collar
v_mod = 100 .* (deg2rad(w) .* r) .* (-sind(theta) - cosd(theta) .* tand(beta));
end