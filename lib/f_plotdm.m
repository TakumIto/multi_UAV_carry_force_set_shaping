function f_plotdm(d, name)

hf=figure;

%subplot(221);
camproj orthographic
h=ezsurf(d, [-pi/2,pi/2,-pi/2,pi/2]);
shading interp
set(h, 'EdgeColor', 'none');
axis vis3d
title('')
x=get(h,'Xdata');
y=get(h,'Ydata');
z=get(h,'Zdata');
hold on
for i = 1:3:length(y)
    Y1(1,1:length(y)) = y(i,1);
    X1=x(i,:);
    Z1=z(i,:);
    line(X1,Y1,Z1,'color','black','LineWidth',1) ;
end
X1=[]; Y1=[]; Z1=[];
for i = 1:3:length(x)
    X1(1:length(x),1) = x(1,i);
    Y1=y(:,i);
    Z1=z(:,i);
    line(X1,Y1,Z1,'color','black','LineWidth',1) ;
end
zlabel(['$' name '$'], 'Interpreter', 'latex');
set(gca,'FontSize', 16)
%{
subplot(223);
camproj orthographic
h=ezsurf(d, [-pi/2,pi/2,-pi/2,pi/2]);
shading interp
set(h, 'EdgeColor', 'none');
view(0,90);
axis square
title('')

vars = symvar(d);

subplot(222);
ezplot(subs(d,vars(1),0),[-pi/2,pi/2]);
grid on
xlabel(['$$\' char(vars(2)) '$$~[rad]'], 'Interpreter', 'latex')
ylabel(['$' name '$'], 'Interpreter', 'latex');
title(['$$\' char(vars(1)) ' = 0$$'], 'Interpreter', 'latex')

subplot(224);
ezplot(subs(d,vars(2),0),[-pi/2,pi/2]);
grid on
xlabel(['$$\' char(vars(1)) '$$~[rad]'], 'Interpreter', 'latex')
ylabel(['$' name '$'], 'Interpreter', 'latex');
title(['$$\' char(vars(2)) ' = 0$$'], 'Interpreter', 'latex')

%}
set(hf,'Position',[100 100 400 400]);
end

