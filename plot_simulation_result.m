function plot_simulation_result(logname, plot_HFS)

load(strcat('simlog/',logname,'.mat'))

reduced_idx = 1:500:end_time*100;
gamma_reduced = gamma_(reduced_idx,:);
z_scale = 1/5;
dist_start = time(find(x(:,1)>=1,1));
dist_end = time(find(x(:,1)>=4,1));

%% 3D
fig_p = figure(1);
set(fig_p, 'Position', [100 100 800 800]);

hold off
p_pos = plot3(x(:,1),x(:,2),x(:,3),'k', 'LineWidth', 2, 'DisplayName', '$[x~y~z]^\top$'); % position
hold on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
p_force = quiver3(x(reduced_idx,1),x(reduced_idx,2),x(reduced_idx,3), ...
                  thrust_ref(reduced_idx,1),thrust_ref(reduced_idx,2),z_scale*thrust_ref(reduced_idx,3), ...
                  0,'Color','red','LineWidth', 2, 'DisplayName', '$f_p^\mathrm{ref}$');
plot3(thrust_ref(:,1)+x(:,1),thrust_ref(:,2)+x(:,2),thrust_ref(:,3)/10+x(:,3),'r--','LineWidth', 2);
grid on; axis equal
view(15,30)
drawnow

if plot_HFS
    M_ftau_tmp = @(gm) M_ftau_fun(gm);
    HFSs = cell(1,length(reduced_idx));
    for j=1:length(reduced_idx)
        i=reduced_idx(j);
        HFS = compute_HFS(M_ftau_tmp,f_max,gamma_(i,:));
        HFSs{j} = HFS;
    
        try
        [k3,~] = convhull(HFS(:,1),HFS(:,2),HFS(:,3),'Simplify',true);
        p_hfs = trisurf(k3,HFS(:,1)+x(i,1),HFS(:,2)+x(i,2),z_scale*HFS(:,3)+x(i,3), ...
            'FaceColor','blue','FaceAlpha',0.1,'EdgeAlpha',0);
        catch
            disp('hoge')
        end
        drawnow
    end
end

xlabel('$x$ [m], $f_{p,x}^\mathrm{ref}$ [N]', 'FontSize', 24, 'Interpreter', 'latex'); 
ylabel('$y$ [m], $f_{p,y}^\mathrm{ref}$ [N]', 'FontSize', 24, 'Interpreter', 'latex');
zlabel('$z$ [m], $f_{p,z}^\mathrm{ref}/10$ [N]', 'FontSize', 24, 'Interpreter', 'latex');
if plot_HFS
    leg = legend([p_pos,p_force,p_hfs],'Location', 'northwest','interpreter', 'latex','NumColumns',1);
    leg.String{3}='HFS';
else
    leg = legend([p_pos,p_force],'Location', 'northwest','interpreter', 'latex','NumColumns',1);
end
leg.FontSize = 22;
ax = gca; ax.LineWidth = 1.2;

%% 2D
fig_p = figure(2);
set(fig_p, 'Position', [100 100 700 700]);

hold off
p_pos = plot(x(:,1),x(:,2),'k', 'LineWidth', 2, 'DisplayName', '$[x~y]^\top$'); % position
hold on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
p_force = quiver(x(reduced_idx,1),x(reduced_idx,2), ...
                 thrust_ref(reduced_idx,1),thrust_ref(reduced_idx,2), ...
                 0,'Color','red','LineWidth', 2, 'DisplayName', '$f_p^\mathrm{ref}$');
plot(thrust_ref(:,1)+x(:,1),thrust_ref(:,2)+x(:,2),'r--','LineWidth', 2);
xlabel('$x$ [m],\hspace{1em} $f_{p,x}^\mathrm{ref}$ [N]', 'FontSize', 24, 'Interpreter', 'latex'); 
ylabel('$y$ [m],\hspace{1em} $f_{p,y}^\mathrm{ref}$ [N]', 'FontSize', 24, 'Interpreter', 'latex');

dist = patch([1 4 4 1], [-5 -5 10 10], 'k', 'FaceAlpha', 0.1, 'EdgeAlpha', 0, 'DisplayName', 'Wind');
[x_wind, y_wind] = meshgrid(1.5:1:3.5,2:2:6);
quiver(x_wind,y_wind,zeros(size(x_wind)),0.5*ones(size(x_wind)),0,"k",'LineWidth',1)
if plot_HFS
    for j=1:length(reduced_idx)
        i=reduced_idx(j);
        HFS = HFSs{j};
        fz_tmp = thrust_ref(i,3);
        fz_face = [10 10 -10 -10; 10 -10 -10 10; fz_tmp*ones(1,4)]';
        try
            S=(intersectionHull('vert',HFS,'vert',fz_face));
            k=convhull(S.vert(:,1),S.vert(:,2));
            p_hfs = patch(S.vert(k,1)+x(i,1),S.vert(k,2)+x(i,2), ...
                'blue','FaceAlpha',0.1,'EdgeAlpha',0.3);
        catch
        end
        drawnow
    end
end

xlim([-2 8])
ylim([-3 8])

if plot_HFS
    leg = legend([p_pos,p_force,p_hfs,dist],'Location', 'northwest','interpreter', 'latex','NumColumns',1);
    leg.String{3}='HFS';
else
    leg = legend([p_pos,p_force,dist],'Location', 'northwest','interpreter', 'latex','NumColumns',1);
end
leg.FontSize = 22;
grid on; axis equal
ax = gca; ax.LineWidth = 1.2;

%% plot position and attitude
fig_p = figure(3);
set(fig_p, 'Position', [705.0000  491.6667  500  280]);
set(gca, 'Position', [0.1800    0.2737    0.7950    0.7213])

x_ref_x = interp1(x_ref.time, x_ref.Data(:,1), time);
x_ref_y = interp1(x_ref.time, x_ref.Data(:,2), time);
x_ref_z = interp1(x_ref.time, x_ref.Data(:,3), time);

hold off
px = plot(time,x(:,1)-x_ref_x,'LineWidth', 2, 'DisplayName', '$x\!-\!x^\mathrm{ref}$'); % position
hold on
py = plot(time,x(:,2)-x_ref_y,'LineWidth', 2, 'DisplayName', '$y\!-\!y^\mathrm{ref}$');
pz = plot(time,x(:,3)-x_ref_z,'LineWidth', 2, 'DisplayName', '$z\!-\!z^\mathrm{ref}$');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
xlim([0 end_time]);
ylim([-0.3 0.3]);
xlabel('$\mathrm{Time\,[sec]}$', 'FontSize', 24, 'Interpreter', 'latex'); 
ylabel('$\mathrm{Position\,[m]}$', 'FontSize', 24, 'Interpreter', 'latex');
leg = legend([px,py,pz],'Location', 'northwest','interpreter', 'latex','NumColumns',3);
leg.FontSize = 20;
grid on
ax = gca; ax.LineWidth = 1.2;

%% 
fig_o = figure(4);
set(fig_o, 'Position', [705.0000  491.6667  500  280]);
set(gca, 'Position', [0.1800    0.2737    0.7950    0.7213])

hold off
pphi = plot(time,(180/pi)*x(:,7),'LineWidth', 2, 'DisplayName', '$\phi$'); % position
hold on
ptheta = plot(time,(180/pi)*x(:,8),'LineWidth', 2, 'DisplayName', '$\theta$');
ppsi = plot(time,(180/pi)*x(:,9),'LineWidth', 2, 'DisplayName', '$\psi$');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
xlim([0 end_time]);
ylim([-1.0 1.0]);
yticks([-0.5 0 0.5])
xlabel('$\mathrm{Time\,[sec]}$', 'FontSize', 24, 'Interpreter', 'latex'); 
ylabel('$\mathrm{Attitude\,[deg]}$', 'FontSize', 24, 'Interpreter', 'latex');
leg = legend([pphi,ptheta,ppsi],'Location', 'northeast','interpreter', 'latex','NumColumns',3);
leg.FontSize = 22;
grid on
ax = gca; ax.LineWidth = 1.2;

%% plot rotor inputs
fig_r = figure(5);
set(fig_r, 'Position', [705.0000  491.6667  500  280]);
set(gca, 'Position', [0.1800    0.2737    0.7950    0.7213])

hold off
plot(time, rotor_sat_input, 'k','LineWidth', 1)
hold on
plot([0 end_time], [0, 0], '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
plot([0 end_time], [f_max, f_max], '--', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
xlim([0 end_time]);
ylim([-0.2 4.2]);
xlabel('$\mathrm{Time\,[sec]}$', 'FontSize', 24, 'Interpreter', 'latex'); 
ylabel('$\mathrm{Rotor\,thrust\,[N]}$', 'FontSize', 22, 'Interpreter', 'latex');
grid on
ax = gca; ax.LineWidth = 1.2;

%% plot angles
fig_a = figure(6);
set(fig_a, 'Position', [705.0000  491.6667  500  280]);
set(gca, 'Position', [0.1900    0.2737    0.7950    0.7213])

hold off
g1 = plot(time, (180/pi)*gamma_(:,2),'LineWidth', 2, 'DisplayName', '$\gamma_1$');
hold on
g2 = plot(time, (180/pi)*gamma_(:,3),'LineWidth', 2, 'DisplayName', '$\gamma_2$');
g3 = plot(time, (180/pi)*gamma_(:,4),'LineWidth', 2, 'DisplayName', '$\gamma_3$');
g4 = plot(time, (180/pi)*gamma_(:,1),'LineWidth', 2, 'DisplayName', '$\gamma_4$');

dist = patch([dist_start dist_end dist_end dist_start], [-50 -50 50 50], 'k', 'FaceAlpha', 0.1, 'EdgeAlpha', 0, 'DisplayName', 'Wind');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
xlim([0 end_time]);
ylim([-24 2]);
xlabel('$\mathrm{Time\,[sec]}$', 'FontSize', 24, 'Interpreter', 'latex'); 
ylabel('$\mathrm{Tilt~angle\,[deg]}$', 'FontSize', 24, 'Interpreter', 'latex');
leg = legend([g1,g2,g3,g4],'Location', 'northeast','interpreter', 'latex','NumColumns',4);
leg.FontSize = 22;
grid on
ax = gca; ax.LineWidth = 1.2;

%% plot x-y force inputs
fig_a = figure(7);
set(fig_a, 'Position', [705.0000  491.6667  500  280]);
set(gca, 'Position', [0.1800    0.2737    0.7950    0.7213])

hold off
set(gca,'ColorOrderIndex',2)
pnomx = plot(Fp_nom.Time, Fp_nom.Data(:,1),'r--','LineWidth', 2, 'DisplayName', '$f_{p,x}^\mathrm{nom}$');
hold on
pnomy = plot(Fp_nom.Time, Fp_nom.Data(:,2),'b--','LineWidth', 2, 'DisplayName', '$f_{p,y}^\mathrm{nom}$');
px = plot(time, thrust_ref(:,1),'r','LineWidth', 2, 'DisplayName', '$f_{p,x}^\mathrm{ref}$');
py = plot(time, thrust_ref(:,2),'b','LineWidth', 2, 'DisplayName', '$f_{p,y}^\mathrm{ref}$');

dist = patch([dist_start dist_end dist_end dist_start], [-5 -5 10 10], 'k', 'FaceAlpha', 0.1, 'EdgeAlpha', 0, 'DisplayName', 'Wind');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);
xlim([0 end_time]);
ylim([-0.7 1.2]);
yticks([-1, 0, 1, 2])
xlabel('$\mathrm{Time\,[sec]}$', 'FontSize', 24, 'Interpreter', 'latex'); 
ylabel('$\mathrm{Force\,[N]}$', 'FontSize', 24, 'Interpreter', 'latex');
leg = legend([pnomx, pnomy, px, py],'Location', 'northwest','interpreter', 'latex','NumColumns',4);
leg.FontSize = 20;
grid on
ax = gca; ax.LineWidth = 1.2;

end
