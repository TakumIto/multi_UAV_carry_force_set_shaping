function [outputArg1,outputArg2] = overlook_platform(N, alphas, betas, gammas, scale)
%OVERLOOK_PLATFORM この関数の概要をここに記述
%   詳細説明をここに記述

plot3(0,0,0);hold on;
axis vis3d; grid on; axis equal;% axis off;
view(0,90)
xlabel('$x\mathrm{[m]}$', 'interpreter', 'latex');
ylabel('$y\mathrm{[m]}$', 'interpreter', 'latex');

%% 
n_rotor = 4;
r_param=scale*0.25;
size_p = scale*[0.7,0.7,0.4]; % payload size x-y-z
size_a = scale*[0.2,0.2,0.05]; % agents' size x-y-z
l = scale*1;
radius = scale*0.1; % propeller radius
direct_length = scale*2;

% agents' rotor pos-ori structure 
alpha_param = repmat(pi/4:2*pi/n_rotor:2*pi-0.0001,N,1); % normal rotor instalation (e.g. [0, pi/2, pi, 3pi/2] for simple-quadrotor)
beta_param = zeros(N,n_rotor);
gamma_param = zeros(N,n_rotor);

% agents position
l_param = l*ones(1,N);
% x_param = l_param.*cos(alphas);
% y_param = l_param.*sin(alphas);
x_param = l_param.*sin(alphas);
y_param = -l_param.*cos(alphas);
z_param = zeros(1,N);

%% object shape definition
cube_face = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];

% payload
p_shape = (1/2).*[-size_p(1) -size_p(2) -size_p(3);size_p(1) -size_p(2) -size_p(3);size_p(1)  size_p(2) -size_p(3);-size_p(1)  size_p(2) -size_p(3);...
                  -size_p(1) -size_p(2)  size_p(3);size_p(1) -size_p(2)  size_p(3);size_p(1)  size_p(2)  size_p(3);-size_p(1)  size_p(2)  size_p(3)];
payload = patch('Vertices',p_shape,'Faces',cube_face,'FaceVertexCData',hsv(1),'FaceColor','#696969');

% agent
a_shape = (1/2).*[-size_a(1) -size_a(2) -size_a(3);size_a(1) -size_a(2) -size_a(3);size_a(1)  size_a(2) -size_a(3);-size_a(1)  size_a(2) -size_a(3);...
                   -size_a(1)  -size_a(2)  size_a(3);size_a(1) -size_a(2)  size_a(3);size_a(1)  size_a(2)  size_a(3);-size_a(1)  size_a(2)  size_a(3)];
% agents = cell(1,n_agent);
for i=1:N
    agents(i) = patch('Vertices',a_shape,'Faces',cube_face,'FaceVertexCData',hsv(1),'FaceColor','#696969');
    directions(i) = plot3(0,0,0,'r','LineWidth',4);
end
for i=1:N
    bar(i) = plot3(0,0,0,'k','LineWidth',2);
end
    
% rotor
n_rotor_vert = 20;
for i = 1:n_rotor_vert
    vert_circle(:,i) = [radius*cos(2*pi/n_rotor_vert*i); radius*sin(2*pi/n_rotor_vert*i); 0];
end
% rotors = cell(n_agent,n_rotor);
for i=1:N
    for j=1:n_rotor
        rotors(i,j) = patch(vert_circle(1,:),vert_circle(2,:),vert_circle(3,:),'FaceColor','#E59100');
    end
end


%% plot
% payload 
pos_p = [0,0,0];
ori_p = [0,0,0];
R_p = Rxyz(ori_p(1),ori_p(2),ori_p(3));
vert_p = (R_p*p_shape')' + repmat(pos_p,8,1);

set(payload,'Vertices',vert_p);

% agent i
for i=1:N
    pos_ai = [x_param(i),y_param(i),z_param(i)];
    ori_ai = [gammas(i),betas(i),alphas(i)];

    R_ai = Rzyx(ori_ai(3),ori_ai(2),ori_ai(1));
    vert_ai = (R_ai*a_shape')' + repmat(pos_ai,8,1);
    vert_ai = (R_p*vert_ai')' + repmat(pos_p,8,1);
    set(agents(i),'Vertices',vert_ai);

    agent_direct = [pos_ai*R_p'+pos_p;[0,0,direct_length]*R_ai'*R_p'+(pos_ai*R_p'+pos_p)];
    set(directions(i),'XData',agent_direct(:,1),'YData',agent_direct(:,2),'ZData',agent_direct(:,3))

    agent_arm = [pos_p;pos_ai*R_p'+pos_p];
    set(bar(i),'XData',agent_arm(:,1),'YData',agent_arm(:,2),'ZData',agent_arm(:,3));

    % rotor j 
    for j=1:n_rotor
        ori_rj = [gamma_param(i,j),beta_param(i,j),alpha_param(i,j)];
        pos_rj = [r_param*cos(ori_rj(3)),r_param*sin(ori_rj(3)),0];
        R_r = Rzyx(ori_rj(3),ori_rj(2),ori_rj(1));
        vert_r = (R_r*vert_circle) + repmat(pos_rj',1,n_rotor_vert);
        vert_r = (R_ai*vert_r) + repmat(pos_ai',1,n_rotor_vert);
        vert_r = (R_p*vert_r) + repmat(pos_p',1,n_rotor_vert);
        set(rotors(i,j),'Vertices',vert_r');
    end
end
end

