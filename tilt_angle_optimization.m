clear; close all;

%% define the parameters 

% load the physical parameters
run('params/physical_parameters.m')

% parameters of RFS
f_c = m*[0.5; 0.5; 9.8]; % center
f_cx_lim = [-1 1];       % fcx range
f_cy_lim = [-1 1];       % fcy range
f_cz_lim = [-1 1];       % fcz range

% optimization options
max_iterations = 30;
max_stall_iterations = 500;
swarm_size = 300; % Particle Swarm Optimization

%% formulation
gamma = sym('gamma',[N 1],'real');

M_quad = [0 0 0 0;
          0 0 0 0;
          1 1 1 1;
          r, r, -r, -r;
          -r, r, r, -r;
          kappa -kappa kappa -kappa];

M_ftau = sym(zeros(6,4*N));
M_gamma = sym(zeros(N, 4*N));
for i=1:N
    R_pqi = Rzyx(alpha_param(i),0,gamma(i));
    l_hat = wedge(p_pqi(:,i));
    M_ftau(:,N*(i-1)+1:N*i) = [R_pqi zeros(3);
                               l_hat*R_pqi R_pqi] * diag([1, 1, 1, 0, 1, 1]) * M_quad;
    
    M_gamma(i,N*(i-1)+1:N*i) = [0, 0, 0, 1, 0, 0] * M_quad;
end
M_ftau = simplify(M_ftau);
M_gamma = double(M_gamma);

% substitution
M_ftau_fun = matlabFunction(M_ftau, 'Vars', {gamma'});

%% Vertices of required force set
cuboid = table2array(combinations(f_cx_lim, f_cy_lim, f_cz_lim));
RFS = f_c'+cuboid;

%% Optimization
try 
    close(1)
    close(2)
catch 
end
ax_fval = axes(figure(1));
ax_x = axes(figure(2));

% settings
limits = [0,f_max];

fun = @(x)-inner_point(double(M_ftau_fun(x)),RFS', limits)...
          + 1/(N*f_max^2)*norm(x)^2; % L2 regularization
nvars = N;
initial_param = [-0.25, -0.25, -0.25, -0.25];
lb = -pi/4*ones(1,nvars);
ub = zeros(1,nvars);

outputfun = @(optimValues, state)plotoutput(optimValues, state, ax_fval, ax_x, nvars);
outputfun_SA = @(options,optimvalues,flag)plotoutput_SA(options,optimvalues,flag, ax_fval, ax_x, nvars);

options = optimoptions('particleswarm','SwarmSize',swarm_size,'MaxIterations',max_iterations,'HybridFcn',@fmincon,...
                       'Display','none','OutputFcn', outputfun,'UseParallel',true,'initialPoints',initial_param,...
                       'FunctionTolerance', 1e-8, 'SelfAdjustmentWeight', 0.1, 'SocialAdjustmentWeight', 0.1,...
                       'MaxStallIterations', max_stall_iterations);

% run the optimization    
tic
opt_angles = particleswarm(fun,nvars,lb,ub,options)
toc

% display the result
fprintf('The number of inner points: %d / %d \n',...
    inner_point(double(M_ftau_fun(opt_angles)),RFS', limits), length(RFS)); 

%% Hoverable Force Set
HFS = compute_HFS(M_ftau_fun,f_max,opt_angles);

%% Visualize
xy_enlarge = 5;

figHandle = figure(3);
set(figHandle,'Position',[100 100 800 800]);
set(gcf,'Renderer','painters');
view(-15,10)

hold off
% plot required force sets
[k1,~] = convhull(RFS(:,1),RFS(:,2),RFS(:,3),'Simplify',true);
trisurf(k1,xy_enlarge*RFS(:,1),xy_enlarge*RFS(:,2),RFS(:,3),'FaceColor','blue','FaceAlpha',0.15,'EdgeAlpha',0.15,'DisplayName','RFS');
% plot hoverable force set
hold on
[k3,~] = convhull(HFS(:,1),HFS(:,2),HFS(:,3),'Simplify',true);
trisurf(k3,xy_enlarge*HFS(:,1),xy_enlarge*HFS(:,2),HFS(:,3),'FaceColor','blue','FaceAlpha',0.1,'EdgeAlpha',0.1,'DisplayName','HFS');

X=[0,0;0,0]; Y = [0,0;0,0]; Z = [0,0;0,0]; % dummy data
p1 = surf(X,Y,Z,'FaceColor',[0.5 0.5 1],'FaceAlpha',0.5,'EdgeAlpha',0.5,'DisplayName','RFS');
p2 = surf(X,Y,Z,'FaceColor',[0.8 0.8 1],'FaceAlpha',0.2,'EdgeAlpha',0.2,'DisplayName','HFS');

% plot platform
visualize_platform(N, alpha_param', zeros(4,1), opt_angles',10)
axis equal

set(figHandle,'color','white'); 
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

xlabel('$f_{px}\mathrm{[N]}$','interpreter', 'latex');
ylabel('$f_{py}\mathrm{[N]}$','interpreter', 'latex');
zlabel('$f_{pz}\mathrm{[N]}$','interpreter', 'latex');

tickvals = [-10, 0, 10];
xticks(tickvals);
yticks(tickvals);
zticks(0:10:70);
xticklabels(string(tickvals/xy_enlarge))
yticklabels(string(tickvals/xy_enlarge))

[~, objh] = legend([p1 p2],'Location', 'north','interpreter', 'latex','NumColumns',2);
objhl = findobj(objh, 'type', 'patch'); % objects of legend of type patch
set(objhl, 'Markersize', 15); % set marker size as desired
view(20,20)

% plot overlooks
figHandle2 = figure(4);
set(figHandle2,'Position',[100 100 800 800]);
hold off
overlook_platform(N, alpha_param', zeros(4,1), opt_angles',10)
hold on
box on

trisurf(k1,xy_enlarge*RFS(:,1),xy_enlarge*RFS(:,2),RFS(:,3),'FaceColor','blue','FaceAlpha',0.15,'EdgeAlpha',0.15,'DisplayName','RFS');
trisurf(k3,xy_enlarge*HFS(:,1),xy_enlarge*HFS(:,2),HFS(:,3),'FaceColor','blue','FaceAlpha',0.1,'EdgeAlpha',0.1,'DisplayName','HFS');

set(figHandle2,'color','white'); 
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

xlabel('$f_{px}\mathrm{[N]}$','interpreter', 'latex');
ylabel('$f_{py}\mathrm{[N]}$','interpreter', 'latex');

tickvals = [-10, 0, 10];
xticks(tickvals);
yticks(tickvals);
xticklabels(string(tickvals/xy_enlarge))
yticklabels(string(tickvals/xy_enlarge))

%% functions 
function n_in = inner_point(M_ftau, RFS, limits)
    n = size(M_ftau,2); % input size (= 4*N)
    
    f = [zeros(1,n), 1];
    A = [-eye(n), limits(1)*ones(n,1);
         eye(n), -limits(2)*ones(n,1);
         zeros(1,n), -1];
    b = zeros(2*n+1,1);
    Aeq = [M_ftau, zeros(6,1)];
    options = optimoptions('linprog','Algorithm','interior-point','Display','none');
    
    n_in = 0;
    parfor i=1:length(RFS) % parallel for
        Freq = RFS(:,i);
        [~, fval] = linprog(f,A,b,Aeq,[Freq;0;0;0],[],[],options);
        if fval <= 1
            n_in = n_in + 1; % step function
        end
    end
end


function stop = plotoutput(optimValues, state, ax_fval, ax_x, nvars)
    bestfval = optimValues.bestfval;
    bestx = optimValues.bestx;
    
    if isempty(ax_fval.Children)
        plot(ax_fval, [bestfval;bestfval])
    else 
        fval_history = ax_fval.Children.YData;
        plot(ax_fval, [fval_history bestfval])
    end

    if isempty(ax_x.Children)
        plot(ax_x, [bestx;bestx])
    else
        for i=1:nvars
            ax_x.Children(i).YData = [ax_x.Children(i).YData, bestx(i)];
        end
    end
    drawnow;
    stop=0;
end

function [stop,options,optchanged] = plotoutput_SA(options,optimvalues,flag, ax_fval, ax_x, nvars)

    stop = false;
    optchanged = false;

    bestfval = optimvalues.bestfval;
    bestx = optimvalues.bestx;
    
    if isempty(ax_fval.Children)
        plot(ax_fval, [bestfval;bestfval])
    else 
        fval_history = ax_fval.Children.YData;
        plot(ax_fval, [fval_history bestfval])
    end

    if isempty(ax_x.Children)
        plot(ax_x, [bestx;bestx])
    else
        for i=1:nvars
            ax_x.Children(i).YData = [ax_x.Children(i).YData, bestx(i)];
        end
    end
    drawnow;
end
