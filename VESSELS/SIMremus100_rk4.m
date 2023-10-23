% SIMremus100_original_rk4.m Vehicle Simulator
clear;
clc;
close all;
% clear all variables
fprintf ('\n\nREMUS DYNAMICS SIMULATOR WITH RK4 INTEGRATION\n');
%disp(sprintf (' NOTE: Model using %s REMUS dimensions. \n\n', vehicle));
%
% Output flags
%
show_step    = 1; 
show_speed   = 0;
show_pos     = 0;
run_savedata = 0; 
run_plots    = 0;
choose_int   = 0;
%choose_setup = 1 ;

% create .mat files
d = clock ; yy = d(1) ; mo = d(2) ; dd = d(3) ; hh = d(4) ; 
mm = d(5) ; ss = d(6) ; 
date_string = datestr (datenum(yy ,mo, dd) ,1);
time_string = datestr(datenum(yy,mo,dd,hh,mm,ss) ,13) ;

% EXPERIMENTAL/ASSIGNED VALUES: initial conditions, input vector
% ------------------------------------------------------------------------------
% loading model inputs, generated in SIM_SETUP. M
% load input_vector % data from FIN_INPUTS. M on mission files
% load time_step
% load initial_state % data from INITIAL_CONDITIONS.M on above

ui = zeros(3,1000);
ui(1,:) = deg2rad(-0.11); % delta_r
ui(2,:) = deg2rad(-2.15); % delta_s
ui(3,:) = 1418; % RPM
time_step = 0.01;
x = [2.0575 0.001 -0.021 0 0 0 0 0 0 deg2rad(-2.42) deg2rad(-0.583) 0]'; % Add trim state
XX = [];

pitch_max = 90 ;

%RUN MODEL
% ------------------------------------------------------------------------------
% Initialize number of steps and storage matrix
n_steps = 2;
% size(ui, 2)-1;
%output_table = zeros(n_steps, size(x,1) + size(ui,1)+7);
fprintf('\nSimulator running...\n');

% MAIN PROGRAM
for i = 1:n_steps
    
    XX = [XX x];
    
    % Print current step for error checking
    if show_step == 1
        if ~rem(i*10 ,n_steps)
             fprintf ('Steps Completed : %02d : %02d%% \n',n_steps,i/n_steps*100);
        end
    end

    % Calculate forces, accelerations
    [xdot, U, Forces, Moments] = remus100_withM(x,ui(:,i)');    

    %% RUNGE-KUTTA APPROXIMATION to calculate new states
    %% NOTE: ideally, should be approximating ui values for k2,k3
    %% ie (ui(: ,i)+ui(: ,i+1))/2
    k1_vec = xdot;
    k2_vec = remus100_withM(x+( 0.5 .* time_step.*k1_vec),((ui(:,i)+ui(:,i+1))./2)');
    k3_vec = remus100_withM(x+( 0.5 .*time_step.*k2_vec),((ui(:,i)+ui(:,i+1))./2)');
    k4_vec = remus100_withM(x+(time_step.*k3_vec) , ui(: ,i+1)') ;
    x = x + time_step/6.*(k1_vec+2.*k2_vec+2.*k3_vec+k4_vec) ;
end


%Plot Output
for i = 1:size(x,1)
    figure;
    plot([1:n_steps],XX(i,:),'*');
    title(sprintf('%d',i));
end

fprintf('Simulation Complete.\n\n');

return;