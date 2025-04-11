clear; close all;

%% physical parameter definition
run('variable_system.m')

%% simulation parameters
logname = 'simulation';

tau_rotor = 0.01; % time constant of rotor
F_dist = [0, 0.5, 0, 0, 0, 0];

% generate LQI gainmatrix
Q=diag([1/10;1/10;400;   
        5;5;10;  
        100;100;100;  
        5;5;10;  
        1;1;1/2;  
        1/10;1/10;1/10]); %state: x y z dx dy dz gamma beta alpha x/s y/s z/s gamma/s beta/s alpha/s
R=diag([100;100;100; 100;100;100]); % input:[F;M]

run("generate_lqr_gainmatrix.m");

% initial states
load(strcat('x_nom.mat'))
load(strcat('Fp_nom.mat'))

x0 = x_ref.Data(1,:)';

% load optimized angle data
load('optimal_tilt_angles.mat')
opt_angle_table_mat = zeros([4,size(opt_angle_table_all{1})]);
gamma0 = zeros(4,1);
for i=1:4
    opt_angle_table_mat(i,:,:) = smoothdata2(opt_angle_table_all{i});
    opt_angle_tmp = reshape(opt_angle_table_mat(i,:,:),size(AX_all));
end

%% simulation
out = sim('multi_UAV_simulation_model.slx');
save(strcat('simlog/',logname,'.mat'))

%% plot result
plot_simulation_result(logname, false)
