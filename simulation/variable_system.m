%% path management
physical_parameters

%% general dynamics
syms t real
syms m_sym positive
J_sym = sym('J',[3 3],'real');

%state
syms x y z phi theta psi real
syms dx dy dz dphi dtheta dpsi real
syms ddx ddy ddz ddphi ddtheta ddpsi real

p_oc=[x;y;z];
R_oc=rot_z(psi)*rot_y(theta)*rot_x(phi);
r_oc=[phi;theta;psi];

% forces of agents to object wrench
F_p = sym('F_p',[6 1],'real');

%% M all
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

matlabFunction(M_ftau,'File','./lib/M_ftau_fun','Vars',{gamma'});

%% EQM
F_c = F_p + [R_oc.'*[0; 0; -m*g]; 0; 0; 0]; % object wrench
f_c = F_c(1:3);
tau_c = F_c(4:6);

v_c = simplify(R_oc.'*tderiv(p_oc, t));
dv_c = simplify(tderiv(v_c, t));

omegaw_c = simplify(R_oc.'*tderiv(R_oc,t));
omega_c = vee(omegaw_c);
domega_c = simplify(tderiv(omega_c, t));

M = [m*eye(3) zeros(3);
     zeros(3) J];
C = [m*omegaw_c*v_c;
     omegaw_c*J*omega_c];

EQM = (simplify(M * [dv_c; domega_c] + C) == F_c);
State_EQM = cell2sym(struct2cell(solve(EQM,[ddx;ddy;ddz;ddphi;ddtheta;ddpsi])));

%% calcurate fx, gx
tmpfx = sym(zeros(length(State_EQM),1));
tmpgx = sym(zeros(length(State_EQM),length(F_p)));

clear tmpc tmpT
for i=1:length(State_EQM)
	[tmpc, tmpT]=coeffs(State_EQM(i),F_p);
	for j=1:length(tmpT)
        if tmpT(j)==1
			tmpfx(i,1)=tmpc(j);
        else
            for k=1:length(F_p)    
                if tmpT(j)==F_p(k)
                    tmpgx(i,k)=tmpc(j);
                end
            end
        end
	end
end
clear tmpc tmpT

fx=[dx;dy;dz;tmpfx(1:3);dphi;dtheta;dpsi;tmpfx(4:6)];
gx=[zeros(3,6);tmpgx(1:3,:);zeros(3,6);tmpgx(4:6,:)];

%% save
statex=[x;y;z;dx;dy;dz;phi;theta;psi;dphi;dtheta;dpsi];
statexd=zeros(12,1);

fx_phys = subs(fx,m_sym,m);
fx_phys = subs(fx_phys,J_sym,J);
gx_phys = subs(gx,m_sym,m);
gx_phys = subs(gx_phys,J_sym,J);
matlabFunction(fx_phys,'File','./lib/fx_fun','Vars',{statex,gamma});
matlabFunction(gx_phys,'File','./lib/gx_fun','Vars',{statex(7:9)});
