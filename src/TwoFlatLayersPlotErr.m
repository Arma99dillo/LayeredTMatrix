function TwoFlatLayersPlotErr(Tmat,Nlayer,b,a_d,a_u)

%plot solution and error against exact solution for the two flat layers problem

%get parameters
param=Tmat{1}.param; param.d1=cos(param.theta); param.d2=sin(param.theta); 
param.gamma=sqrt(param.epsilon(2)-param.d1^2); param.L1=(Tmat{1}.Hp+Tmat{1}.Hm)/2;

%compute exact solution given the parameters
P = [exp(-1i*param.K*param.L1*param.d2), -exp(-1i*param.K*param.L1*param.gamma), -exp(1i*param.K*param.L1*param.gamma), 0;
    -param.d2*exp(-1i*param.K*param.L1*param.d2), param.gamma*exp(-1i*param.K*param.L1*param.gamma), -param.gamma*exp(1i*param.K*param.L1*param.gamma), 0;
    0, exp(1i*param.K*param.L1*param.gamma), exp(-1i*param.K*param.L1*param.gamma), -exp(-1i*param.K*param.L1*param.d2);
    0, -param.gamma*exp(1i*param.K*param.L1*param.gamma), param.gamma*exp(-1i*param.K*param.L1*param.gamma), -param.d2*exp(-1i*param.K*param.L1*param.d2)];
l = [-exp(1i*param.K*param.L1*param.d2); -param.d2*exp(1i*param.K*param.L1*param.d2); 0; 0];

V = P\l; param.R=V(1); param.T1=V(2); param.T2=V(3); param.T3=V(4);

uex = {@(x1,x2) exp(1i*param.K.*(x1.*param.d1+x2.*param.d2)) + ...
        param.R.*exp(1i*param.K.*(x1.*param.d1-x2.*param.d2)), ...
        @(x1,x2) param.T1.*exp(1i*param.K.*(x1.*param.d1-x2.*param.gamma)) + ...
        param.T2.*exp(1i*param.K.*(x1.*param.d1+x2.*param.gamma)), ...
        @(x1,x2) param.T3.*exp(1i*param.K.*(x1.*param.d1+x2.*param.d2))};

%plot solution and error
figure()
disp('Computing and plotting the T-matrix solution between the layers')
tic()

M=param.M; L=param.L; alpha_n=Tmat{1}.alpha_n; Ndim=2*M+1;

%functions Phi_n to evaluate the T-matrix solution
eval_up = @(x1,x2,alpha,beta,Hp) exp(1i.*beta.*(x2-Hp)).*exp(1i.*alpha.*x1);
eval_down = @(x1,x2,alpha,beta,Hm) exp(-1i.*beta.*(x2-Hm)).*exp(1i.*alpha.*x1);

%plot solution and error between the layers
r=1;
for j=1:Nlayer-1

    %setup grid in the middle
    npoints=200;
    g1=linspace(0,L,npoints); g2=linspace(Tmat{j+1}.Hp,Tmat{j}.Hm,npoints);
    [x1,x2]=meshgrid(g1,g2);
    
    %get vector for upward and downward expansion
    b_d=b(r*Ndim+1:(r+1)*Ndim); b_u=b((r+1)*Ndim+1:(r+2)*Ndim); 

    %compute solution using the T-matrix expansion
    sol=zeros(npoints,npoints); v=1;
    for n=-M:M
        sol=sol+b_u(v).*eval_up(x1,x2,alpha_n(v),Tmat{j+1}.beta_n_p(v),Tmat{j+1}.Hp);
        sol=sol+b_d(v).*eval_down(x1,x2,alpha_n(v),Tmat{j}.beta_n_m(v),Tmat{j}.Hm);
        v=v+1;
    end
    uex_mid=uex{2}(x1,x2);

    %plot solution and error against exact solution
    subplot(1,3,1)
    surf(x1,x2,real(sol)); shading interp;
    view(2);
    hold on

    subplot(1,3,2)
    surf(x1,x2,abs(sol));  shading interp;
    view(2);
    hold on

    subplot(1,3,3)
    surf(x1,x2,log10(abs(sol-uex_mid))); shading interp;
    view(2);
    hold on
    
    if j == 1
        minColorLimit1 = min(min(real(sol)));  % determine colorbar limits from data
        maxColorLimit1 = max(max(abs(sol)));
        minColorLimit2 = min(min((log10(abs(sol-uex_mid)))));
        maxColorLimit2 = max(max(log10(abs(sol-uex_mid))));
    else
        minColorLimit1 = min(min(min(real(sol))),minColorLimit1);  % determine colorbar limits from data
        maxColorLimit1 = max(max(max(abs(sol))),maxColorLimit1);
        minColorLimit2 = min(min(min((log10(abs(sol-uex_mid))))),minColorLimit2);
        maxColorLimit2 = max(max(max(log10(abs(sol-uex_mid)))),maxColorLimit2);
    end

    r=r+2;
end

toc()

%then plot solution and error inside layers

%first layer: get exact solution and translate mesh
disp('Plot inside layer 1')
uex_temp = {@(x) uex{1}(x(1),x(2)), @(x) uex{2}(x(1),x(2))};
mesh_temp=Tmat{1}.mesh;
mesh_temp.p(:,2) = mesh_temp.p(:,2)+(Tmat{1}.Hp-Tmat{1}.param.H); %translate mesh

%get T-matrix expansion vectors
b_u_p=b(1:Ndim); b_d_p=a_d; b_d_m=b(Ndim+1:2*Ndim); b_u_m=b(2*Ndim+1:3*Ndim);
deltaH_p=0; deltaH_m=Tmat{1}.Hm-Tmat{2}.Hp;

%solve impedance boundary value problem
A = MatrixTDGImp(mesh_temp,Tmat{1}.param);
rhs = rhsTDGImp(mesh_temp,Tmat{1},b_u_p,b_d_p,b_u_m,b_d_m,deltaH_p,deltaH_m);
u = A\rhs;

%plot solution and error against exact solution
[minColorLimit1,maxColorLimit1,minColorLimit2,maxColorLimit2]=...
    PlotErrorLog(mesh_temp,Tmat{1}.param,u,uex_temp,minColorLimit1,maxColorLimit1,minColorLimit2,maxColorLimit2);
hold on

%last layer: get exact solution and translate mesh
disp(['Plot inside layer ', num2str(Nlayer)])
uex_temp = {@(x) uex{2}(x(1),x(2)), @(x) uex{3}(x(1),x(2))};
mesh_temp=Tmat{Nlayer}.mesh;
mesh_temp.p(:,2) = mesh_temp.p(:,2)+(Tmat{Nlayer}.Hp-Tmat{Nlayer}.param.H); %translate mesh

%get T-matrix expansion vectors
b_d_p=b(2*Nlayer*Ndim-3*Ndim+1:2*Nlayer*Ndim-2*Ndim); b_u_p=b(2*Nlayer*Ndim-2*Ndim+1:2*Nlayer*Ndim-Ndim);
b_d_m=b(2*Nlayer*Ndim-Ndim+1:end); b_u_m=a_u;
deltaH_p=Tmat{Nlayer-1}.Hm-Tmat{Nlayer}.Hp; deltaH_m=0;

%solve impedance boundary value problem
A = MatrixTDGImp(mesh_temp,Tmat{Nlayer}.param);
rhs = rhsTDGImp(mesh_temp,Tmat{Nlayer},b_u_p,b_d_p,b_u_m,b_d_m,deltaH_p,deltaH_m);
u = A\rhs;

%plot solution and error against exact solution
PlotErrorLog(mesh_temp,Tmat{Nlayer}.param,u,uex_temp,minColorLimit1,maxColorLimit1,minColorLimit2,maxColorLimit2);
hold on

%plot the layers limits and scatterers profile
d=(Tmat{1}.Hp+Tmat{1}.Hm)/2;
subplot(1,3,1)
xlim([0,L]); ylim([Tmat{2}.Hm,Tmat{1}.Hp]);
%plot layers
plot3([0,L],[d,d],10*ones(1,2),'r-',LineWidth=1.8)
plot3([0,L],[Tmat{1}.Hm,Tmat{1}.Hm],10*ones(1,2),'r--',LineWidth=1.8)
plot3([0,L],[Tmat{1}.Hp,Tmat{1}.Hp],10*ones(1,2),'r--',LineWidth=1.8)
plot3([0,L],[-d,-d],10*ones(1,2),'r-',LineWidth=1.8)
plot3([0,L],[Tmat{2}.Hp,Tmat{2}.Hp],10*ones(1,2),'r--',LineWidth=1.8)
plot3([0,L],[Tmat{2}.Hm,Tmat{2}.Hm],10*ones(1,2),'r--',LineWidth=1.8)
subplot(1,3,2)
xlim([0,Tmat{1}.param.L]); ylim([Tmat{2}.Hm,Tmat{1}.Hp]);
%plot layers
plot3([0,L],[d,d],10*ones(1,2),'r-',LineWidth=1.8)
plot3([0,L],[Tmat{1}.Hm,Tmat{1}.Hm],10*ones(1,2),'r--',LineWidth=1.8)
plot3([0,L],[Tmat{1}.Hp,Tmat{1}.Hp],10*ones(1,2),'r--',LineWidth=1.8)
plot3([0,L],[-d,-d],10*ones(1,2),'r-',LineWidth=1.8)
plot3([0,L],[Tmat{2}.Hp,Tmat{2}.Hp],10*ones(1,2),'r--',LineWidth=1.8)
plot3([0,L],[Tmat{2}.Hm,Tmat{2}.Hm],10*ones(1,2),'r--',LineWidth=1.8)

subplot(1,3,3)
xlim([0,Tmat{1}.param.L]); ylim([Tmat{2}.Hm,Tmat{1}.Hp]);
%plot layers
plot3([0,L],[d,d],10*ones(1,2),'r-',LineWidth=1.8)
plot3([0,L],[Tmat{1}.Hm,Tmat{1}.Hm],10*ones(1,2),'r--',LineWidth=1.8)
plot3([0,L],[Tmat{1}.Hp,Tmat{1}.Hp],10*ones(1,2),'r--',LineWidth=1.8)
plot3([0,L],[-d,-d],10*ones(1,2),'r-',LineWidth=1.8)
plot3([0,L],[Tmat{2}.Hp,Tmat{2}.Hp],10*ones(1,2),'r--',LineWidth=1.8)
plot3([0,L],[Tmat{2}.Hm,Tmat{2}.Hm],10*ones(1,2),'r--',LineWidth=1.8)



end