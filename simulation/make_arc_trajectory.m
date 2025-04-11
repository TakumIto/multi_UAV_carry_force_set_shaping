physical_parameters
PI2 = 2*pi;

%% define trajectory
ts = 1/120; % 100Hz

% takeoff
t_takeoff = 10; % time to takeoff
height = 2; % height of arc trajectory
t_stay1 = 5; % time until starting accel

% accelerate x+ direction
t_accel = 5; % accelerating time
v_arc = pi/4; % final velocity

% curve with clothoid path
t_clotho = 5;
arc_start = 0;
arc_end = pi/2;

% decelerate until v=0
t_decel = 5;

% landing
t_stay2 = 5; % time between arc and 
t_landing = 5;

% generate trajectory
t_end = 0;
time = [];
acc_ref = [];
vel_ref = [];
pos_ref = [];

% takeoff
time_tmp = (0:ts:t_takeoff)';
time = [time; t_end+time_tmp];
t_end = t_end+t_takeoff;

acc_ref = [acc_ref; 
           zeros(length(time_tmp),2), height/PI2*((PI2/t_takeoff)^2*sin(PI2/t_takeoff*time_tmp))];
vel_ref = [vel_ref;
           zeros(length(time_tmp),2), height/PI2*(-PI2/t_takeoff*cos(PI2/t_takeoff*time_tmp)+PI2/t_takeoff)];
pos_ref = [pos_ref;
           zeros(length(time_tmp),2), height/PI2*(-sin(PI2/t_takeoff*time_tmp)+PI2/t_takeoff*time_tmp)];

% stay
time_tmp = (0:ts:t_stay1)';
time = [time; t_end+time_tmp];
t_end = t_end+t_stay1;

acc_ref = [acc_ref; 
           zeros(length(time_tmp),3)];
vel_ref = [vel_ref; 
           zeros(length(time_tmp),3)];
pos_end = pos_ref(end,:);
pos_ref = [pos_ref; 
           repmat(pos_end,length(time_tmp),1)];

% accelerating
time_tmp = (0:ts:t_accel)';
time = [time; t_end+time_tmp];
t_end = t_end+t_accel;

acc_ref = [acc_ref; 
           v_arc/2*(pi/t_accel*sin(pi/t_accel*time_tmp)), zeros(length(time_tmp),2)];
vel_ref = [vel_ref;
           v_arc/2*(-cos(pi/t_accel*time_tmp)+1), zeros(length(time_tmp),2)];
pos_end = pos_ref(end,:);
pos_ref = [pos_ref;
           pos_end+[v_arc/2*(-t_accel/pi*sin(pi/t_accel*time_tmp)+time_tmp), zeros(length(time_tmp),2)]];

% clotho
time_tmp = (0:ts:t_clotho)';
time = [time; t_end+time_tmp];
t_end = t_end+t_clotho;

theta = linspace(arc_start,arc_end,length(time_tmp))';
acc_ref = [acc_ref; 
           -v_arc*2*(sqrt(pi/4)/t_clotho)^2*time_tmp.*sin((pi/4/t_clotho)^2*time_tmp.^2), v_arc*2*(sqrt(pi/4)/t_clotho)^2*time_tmp.*cos((pi/4/t_clotho)^2*time_tmp.^2), zeros(length(time_tmp),1)];
vel_ref = [vel_ref;
           v_arc*cos((sqrt(pi/4)/t_clotho)^2*time_tmp.^2), v_arc*sin((sqrt(pi/4)/t_clotho)^2*time_tmp.^2), zeros(length(time_tmp),1)];
init_idx = length(pos_ref);
for i=length(pos_ref)+1:length(vel_ref)
    pos_ref(i,:) = pos_ref(i-1,:)+ts*vel_ref(i-1,:);
end

% duplicate trajectory
time = [time;
        t_end+time];
acc_ref = [acc_ref; 
           -flip(acc_ref(:,2)),-flip(acc_ref(:,1)),flip(acc_ref(:,3))];
vel_ref = [vel_ref; 
           flip(vel_ref(:,2)),flip(vel_ref(:,1)),-flip(vel_ref(:,3))];
pos_end = pos_ref(end,:);
pos_ref = [pos_ref;
           pos_end(1)+flip(-pos_ref(:,2)+pos_end(2)), pos_end(2)+flip(-pos_ref(:,1)+pos_end(1)), flip(pos_ref(:,3))];

% remove same time point
for t=2:length(time)
    if time(t-1)==time(t)
        time(t) = NaN;
        acc_ref(t,:) = NaN(1,3);
        vel_ref(t,:) = NaN(1,3);
        pos_ref(t,:) = NaN(1,3);
    end
end
time = rmmissing(time);
acc_ref = rmmissing(acc_ref);
vel_ref = rmmissing(vel_ref);
pos_ref = rmmissing(pos_ref);

%% plot
fig_p = figure('Position',[0 0 1500 400]);
set(fig_p,'color','white'); 
tiledlayout(1,3,'TileSpacing','tight');

nexttile
hold on; grid on;
plot(time,m*acc_ref, 'LineWidth', 2)
xlim([0 end_time]);
ylim([-0.7 1.2])
xticks(0:10:50)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
xlabel('$\mathrm{Time\,[sec]}$', 'FontSize', 24, 'Interpreter', 'latex'); 
ylabel('$\mathrm{Nominal~force\,[N]}$', 'FontSize', 24, 'Interpreter', 'latex');
leg = legend({"$f_{p,x}^\mathrm{nom}$","$f_{p,y}^\mathrm{nom}$","$f_{p,z}^\mathrm{nom}\!\!-\!\!mg$"}, ...
    'Location', 'northeast','interpreter', 'latex');
leg.FontSize = 22;
ax = gca; ax.LineWidth = 1.2;

nexttile
hold on; grid on;
plot(time,vel_ref, 'LineWidth', 2)
xlim([0 end_time]);
xticks(0:10:50)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
xlabel('$\mathrm{Time\,[sec]}$', 'FontSize', 24, 'Interpreter', 'latex'); 
ylabel('$\mathrm{Velocity\,[m/s]}$', 'FontSize', 24, 'Interpreter', 'latex');
leg = legend({"$v_{x}^\mathrm{ref}$","$v_{y}^\mathrm{ref}$","$v_{z}^\mathrm{ref}$"}, ...
    'Location', 'northeast','interpreter', 'latex');
leg.FontSize = 22;
ax = gca; ax.LineWidth = 1.2;

nexttile
hold on; grid on;
plot(time,pos_ref, 'LineWidth', 2)
xlim([0 end_time]);
xticks(0:10:50)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
xlabel('$\mathrm{Time\,[sec]}$', 'FontSize', 24, 'Interpreter', 'latex'); 
ylabel('$\mathrm{Position\,[m]}$', 'FontSize', 24, 'Interpreter', 'latex');
leg = legend({"$x^\mathrm{ref}$","$y^\mathrm{ref}$","$z^\mathrm{ref}$"}, ...
    'Location', 'northwest','interpreter', 'latex');
leg.FontSize = 22;
ax = gca; ax.LineWidth = 1.2;

%% prepare path reference

% accelaration reference
Fp = [m*(acc_ref+[0 0 g]), zeros(length(time), 3)];
x_all = [pos_ref, vel_ref, zeros(length(time), 6)];
Fp_nom = timeseries(Fp, time);
x_ref = timeseries(x_all,time);

% save
save(strcat('data/Fp_nom'),'Fp_nom','-v7.3');
save(strcat('data/x_nom'),'x_ref','-v7.3');
