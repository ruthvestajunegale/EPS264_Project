clc;clear; 
% clearvars -except del
close all

global Te0 Drho step zf Dflux Dflux_der

% Parameters
rhoprofs = table2array(readtable("rhoprofiles_1to5_030926.csv")); 
muprofs = table2array(readtable("muprofiles_1to5_030926.csv")); 
N = size(rhoprofs,1)-1; %Spatial resolution
zf = linspace(0,1,N+1)-0.5;
step = 0.5;
if step == 1 
    R = 0.9; %Density contrast between core and annulus
    M = 100; %Viscosity contrast between core and annulus
    Drho = -.07/R; %Drop in core density at the nucleation depth, z=0
    Dmu = 9; %Increase in core viscosity at the nuclation depth, z=0
    drho = Drho*heaviside(zf);
    dmu = Dmu*heaviside(zf);
elseif step < 1 && step > 0
    sind =2; 
    R = rhoprofs(1,sind); 
    M = muprofs(1,sind); 
    drho = rhoprofs(:,sind)-R; 
    dmu =  muprofs(:,sind)-M; 
    % dfdx = sech(zf).^2;
    Dflux = zeros(length(zf),1);
    Dflux_der = Dflux;
else 
    R = 0.9; %Density contrast between core and annulus
    M = 100; %Viscosity contrast between core and annulus
    Drho = -.07/R; %Drop in core density at the nucleation depth, z=0
    Dmu = 9; %Increase in core viscosity at the nuclation depth, z=0
    lrho = 0.05;
    drho = Drho*(1+tanh(zf/lrho))/2;
    dfdx = sech(zf/lrho).^2/lrho;
    lmu = 0.05;
    dmu = Dmu*(1+tanh(zf/lmu))/2;
    Dflux = zeros(length(zf),1);
    Dflux_der = Dflux;
end

phi = M./(1+dmu); %Modified viscosity contrast
psi = 1-R*drho/(1-R); %Rescaling factor for the driving force, P

utest = (.00:.01:1).^2;
Ftest = utest;
DFtest = utest;
for i = 1:length(utest)
    Ftest(i) = flux_xlinear(zf(2),utest(i),phi(2),psi(2));
    DFtest(i) = flux_der_xlinear(zf(2),utest(i),phi(2),psi(2));
end

[Ffl,ifl] = max(Ftest);
u_max = max(abs(DFtest)); 
delta_z = 1/N;
zc = zf(2:end)-0.5*delta_z;
CFL = 0.1; 
t_fin = 1;
delta_t = CFL*delta_z/u_max;

% Set initial delta value
u0 = .515^2; 
Tel = flux_xlinear(zf(2),u0,phi(2),psi(2));
Teu = flux_xlinear(zf(2),u0,phi(end-1),psi(end-1));
if abs(Drho) < 1e-6
    Te0 = 0;
else
    Te0 = Tel*heaviside(utest(ifl)-u0)+min(Teu,Ffl)*heaviside(u0-utest(ifl));
end

t = 0:delta_t:t_fin;

% Initialization
u = u0'.*ones(length(zc),1);
u_plt = ones(length(zc),10);
u_plt(:,1) = u;
u_new = u;
a_p12 = zeros(length(zc),1);
a_m12 = zeros(length(zc),1);
f_p12 = zeros(length(zc),1);
f_m12 = zeros(length(zc),1);
u_z = zeros(length(zc),1);
dt_check = 0;

% Explicit high resolution TVD

n = 2;
t_plt = .88/2;
plt = 1;
dt_plt = t_plt;
time = 0;
t_fin = 10;

limit = 9e-4;

while time < t_fin
    fprintf('time: %.3f, new dt: %.2e\n',time,delta_t);
    
    % if step < 0.5
    %     uf = interp1(zc,u,zf);
    %     uf(1) = uf(2);
    %     uf(end) = uf(end-1);
    %     intg = dfdx.*uf.*((2-uf.*(2-phi)).^2./...
    %         (1-uf.^2.*(1-phi))-phi+4*log(uf.^.5));
    %     cintl = cumtrapz(zf,intg);
    %     Dflux = Drho*cintl.*(1-uf.^2).^2/8;
    %     Dflux_der = -Drho*cintl.*(1-uf.^2).*uf/2;
    % end
    
    % Compute flux limiter at time step n
    for j = 2:length(zc)-1
        a = (u(j)-u(j-1))/delta_z;
        b = (u(j+1)-u(j))/delta_z;
        u_z(j) = 0.5*(sign(a)+sign(b))*min(abs(a),abs(b));
    end

    for j = 3:length(zc)-2
        z_p = zf(j+1);
        z_m = zf(j);
        
        phi_cp = phi(j+1);
        psi_cp = psi(j+1);
        phi_cm = phi(j);
        psi_cm = psi(j);

        u_R_p12 = (u(j+1)-delta_z/2*u_z(j+1));
        u_L_p12 = (u(j)+delta_z/2*u_z(j));
        u_R_m12 = (u(j)-delta_z/2*u_z(j));
        u_L_m12 = (u(j-1)+delta_z/2*u_z(j-1));

        a_p12(j) = max( abs(flux_der_xlinear(z_p,u_R_p12,phi_cp,psi_cp)) ,abs(flux_der_xlinear(z_p,u_L_p12,phi_cp,psi_cp)) );
        a_m12(j) = max( abs(flux_der_xlinear(z_m,u_R_m12,phi_cm,psi_cm)) ,abs(flux_der_xlinear(z_m,u_L_m12,phi_cm,psi_cm)) );
        
        f_p12(j) = 0.5*((flux_xlinear(z_p,u_R_p12,phi_cp,psi_cp)+flux_xlinear(z_p,u_L_p12,phi_cp,psi_cp))-...
                            a_p12(j)*(u_R_p12-u_L_p12));
        f_m12(j) = 0.5*((flux_xlinear(z_m,u_R_m12,phi_cm,psi_cm)+flux_xlinear(z_m,u_L_m12,phi_cm,psi_cm))-...
                            a_m12(j)*(u_R_m12-u_L_m12));
                                      
        u_new(j) = u(j)-delta_t/delta_z*(f_p12(j)-f_m12(j));
    end
   
    % Boundary conditions

    u_new(1) = u_new(3);
    u_new(2) = u_new(3); 
    u_new(end) = u_new(end-2);
    u_new(end-1) = u_new(end-2);
    
    if plt
        plot(zc,u_new.^.5,'k','linewidth',2);
        ylim([0 1])
        drawnow
    end
    
    % Update time
    time = time+delta_t;
    
    % CFL e new time step
    dt_check = CFL*delta_z/max(max(a_p12),max(a_m12));
    delta_t = dt_check;
    
    % Update central core flux
%     uc0 = interp1(zc,u_new,0);
%     if step > 0.5
%         Te0 = flux_xlinear(0,uc0,phi(N/2+1),psi(N/2+1));
%     end
    
    % Increase time step
    if time > t_plt
        u_plt(:,n) = u_new;
        t_plt = t_plt+dt_plt;
        n = n+1;
    end
    
    disp("Convergence: " + string(sqrt(mean(((u_new(u_new>0)-u(u_new>0))./u_new(u_new>0)).^2))));
    if sqrt(mean(((u_new(u_new>0)-u(u_new>0))./u_new(u_new>0)).^2))<limit
        disp("converged")
        break
    end

    disp("End stability: " + string(abs(u_new(round(9*N/10-2))-u(round(9*N/10-2)))/u_new(round(9*N/10-2))));
    if (abs(u_new(round(9*N/10-2))-u(round(9*N/10-2)))/u_new(round(9*N/10-2)) > limit) 
        disp("end changes big")
        break
    end

    disp("Beg stability: " + string(abs(u_new(round(3+N/10))-u(round(3+N/10)))/u_new(round(3+N/10))));
    if (abs(u_new(round(3+N/10))-u(round(3+N/10)))/u_new(round(3+N/10)) > limit)
        disp("beg changes big")
        break
    end
    
    u = u_new;
end

filename = ['Dr0Dm9D',num2str(u0^.5*100,'%.0f'),'.mat'];
save(filename,'u_plt','dt_plt');

figure 
plot(zf, drho)

function H = heaviside(z)
H = ones(size(z));
H(z<=0) = 0;
end