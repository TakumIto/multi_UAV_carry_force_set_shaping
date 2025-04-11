% the number of agents and rotors
N = 4; % the number of agent(vehicle)
n_rotor = 4; % the rotors of each vehicle
n = N * n_rotor; % the total number of rotor

% phisical params
g=9.8;
m = 2.5;
J = diag([1/20,1/20,1/20]); % Inertial tensor of payload 
J_xi = diag([1/200,1/200,1/200]); % inertial tensor of agent
kappa=0.015; % counter torque constant
f_max=4; % maximum force of rotor
l=0.22; % length between payload and vehicle 
r=0.08; % distance between rotor and CoG of agent

% agents posision
alpha_param = 0:2*pi/N:2*pi-1e-4; % place the vehicles to uniform interval
for i=1:N
    % pxyz(:,i) = rot_z(alpha_param(i))*[l;0;0];
    p_pqi(:,i) = rot_z(alpha_param(i))*[0;-l;0];
end

end_time = 50; % end time of simulation
