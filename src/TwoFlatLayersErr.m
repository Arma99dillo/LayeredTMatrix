function [ErrL2,ErrH1] = TwoFlatLayersErr(Tmat,Nlayer,b,a_d,a_u)

%compute relative error against exact solution for the two flat layers problem

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

uexdx = {@(x1,x2) 1i*param.K*param.d1.*exp(1i*param.K.*(x1.*param.d1+x2.*param.d2)) + ...
        1i*param.K*param.d1.*param.R.*exp(1i*param.K.*(x1.*param.d1-x2.*param.d2)), ...
        @(x1,x2) 1i*param.K*param.d1.*param.T1.*exp(1i*param.K.*(x1.*param.d1-x2.*param.gamma)) + ...
        1i*param.K*param.d1.*param.T2.*exp(1i*param.K.*(x1.*param.d1+x2.*param.gamma)), ...
        @(x1,x2) 1i*param.K*param.d1.*param.T3.*exp(1i*param.K.*(x1.*param.d1+x2.*param.d2))};

uexdy = {@(x1,x2) 1i*param.K*param.d2.*exp(1i*param.K.*(x1.*param.d1+x2.*param.d2)) - ...
    1i*param.K*param.d2.*param.R.*exp(1i*param.K.*(x1.*param.d1-x2.*param.d2)), ...
    @(x1,x2) -1i*param.K*param.gamma.*param.T1.*exp(1i*param.K.*(x1.*param.d1-x2.*param.gamma)) + ...
    1i*param.K*param.gamma.*param.T2.*exp(1i*param.K.*(x1.*param.d1+x2.*param.gamma)), ...
    @(x1,x2) 1i*param.K*param.d2.*param.T3.*exp(1i*param.K.*(x1.*param.d1+x2.*param.d2))};

%compute error against exact solution
disp('Computing error against exact solution')
tic()

M=param.M; L=param.L; alpha_n=Tmat{1}.alpha_n;
Ndim=2*M+1;

%functions Phi_n to evaluate the T-matrix solution
eval_up = @(x1,x2,alpha,beta,Hp) exp(1i.*beta.*(x2-Hp)).*exp(1i.*alpha.*x1);
eval_down = @(x1,x2,alpha,beta,Hm) exp(-1i.*beta.*(x2-Hm)).*exp(1i.*alpha.*x1);

semil2=0; semih1=0; NormExL2=0; NormExH1=0;

%first compute error between the layers
r=1;
for j=1:Nlayer-1
    
    %setup grid between the layers
    npoints=200;
    grid1=linspace(0,L,npoints); grid2=linspace(Tmat{j+1}.Hp,Tmat{j}.Hm,npoints);
    [x1,x2]=meshgrid(grid1,grid2);
    dx = grid1(2) - grid1(1); dy = grid2(2) - grid2(1);
    
    %get vector for upward and downward expansion
    b_d=b(r*Ndim+1:(r+1)*Ndim); b_u=b((r+1)*Ndim+1:(r+2)*Ndim); 
    
    %compute solution using the T-matrix expansion
    sol=zeros(npoints,npoints); g1=zeros(npoints,npoints); g2=zeros(npoints,npoints); v=1;
    for n=-M:M
        sol=sol+b_u(v).*eval_up(x1,x2,alpha_n(v),Tmat{j+1}.beta_n_p(v),Tmat{j+1}.Hp);
        sol=sol+b_d(v).*eval_down(x1,x2,alpha_n(v),Tmat{j}.beta_n_m(v),Tmat{j}.Hm);
        g1=g1+1i.*alpha_n(v).*b_u(v).*eval_up(x1,x2,alpha_n(v),Tmat{j+1}.beta_n_p(v),Tmat{j+1}.Hp);
        g1=g1+1i.*alpha_n(v).*b_d(v).*eval_down(x1,x2,alpha_n(v),Tmat{j}.beta_n_m(v),Tmat{j}.Hm);
        g2=g2+1i.*Tmat{j+1}.beta_n_p(v).*b_u(v).*eval_up(x1,x2,alpha_n(v),Tmat{j+1}.beta_n_p(v),Tmat{j+1}.Hp);
        g2=g2-1i.*Tmat{j}.beta_n_m(v).*b_d(v).*eval_down(x1,x2,alpha_n(v),Tmat{j}.beta_n_m(v),Tmat{j}.Hm);
        v=v+1;
    end
    sol_ex=uex{2}(x1,x2); g1ex=uexdx{2}(x1,x2); g2ex=uexdy{2}(x1,x2);

    %compute error against exact solution
    semil2 = semil2 + sum(abs(sol - sol_ex).^2, 'all') * dx * dy;
    semih1 = semih1 + sum(abs(sol - sol_ex).^2, 'all') * dx * dy + sum(abs(g1 - g1ex).^2, 'all') * dx * dy + sum(abs(g2 - g2ex).^2, 'all') * dx * dy;
    NormExL2 = NormExL2 + sum(abs(sol_ex).^2, 'all') * dx * dy;
    NormExH1 = NormExH1 + sum(abs(sol_ex).^2, 'all') * dx * dy + sum(abs(g1ex).^2, 'all') * dx * dy + sum(abs(g2ex).^2, 'all') * dx * dy;

    r=r+2;
end

%then compute the error inside layers

%first layer: get exact solution and translate mesh
uex_temp = {@(x) uex{1}(x(1),x(2)), @(x) uex{2}(x(1),x(2))};
uexdx_temp = {@(x) uexdx{1}(x(1),x(2)), @(x) uexdx{2}(x(1),x(2))};
uexdy_temp = {@(x) uexdy{1}(x(1),x(2)), @(x) uexdy{2}(x(1),x(2))};

mesh_temp=Tmat{1}.mesh;
mesh_temp.p(:,2) = mesh_temp.p(:,2)+(Tmat{1}.Hp-Tmat{1}.param.H); %translate mesh

%get T-matrix expansion vectors
b_u_p=b(1:Ndim); b_d_p=a_d; b_d_m=b(Ndim+1:2*Ndim); b_u_m=b(2*Ndim+1:3*Ndim);
deltaH_p=0; deltaH_m=Tmat{1}.Hm-Tmat{2}.Hp;

%solve impedance boundary value problem
A = MatrixTDGImp(mesh_temp,Tmat{1}.param);
rhs = rhsTDGImp(mesh_temp,Tmat{1},b_u_p,b_d_p,b_u_m,b_d_m,deltaH_p,deltaH_m);
u = A\rhs;

%compute error against exact solution
[errl2,norml2,errh1,normh1]=InsideError(mesh_temp,Tmat{1}.param,u,uex_temp,uexdx_temp,uexdy_temp);
semil2=semil2+errl2;
NormExL2=NormExL2+norml2;
semih1=semih1+errh1;
NormExH1=NormExH1+normh1;

%last layer: get exact solution and translate mesh
uex_temp = {@(x) uex{2}(x(1),x(2)), @(x) uex{3}(x(1),x(2))};
uexdx_temp = {@(x) uexdx{2}(x(1),x(2)), @(x) uexdx{3}(x(1),x(2))};
uexdy_temp = {@(x) uexdy{2}(x(1),x(2)), @(x) uexdy{3}(x(1),x(2))};

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

%compute error against exact solution
[errl2,norml2,errh1,normh1]=InsideError(mesh_temp,Tmat{Nlayer}.param,u,uex_temp,uexdx_temp,uexdy_temp);
semil2=semil2+errl2;
NormExL2=NormExL2+norml2;
semih1=semih1+errh1;
NormExH1=NormExH1+normh1;

toc()

%return relative errors
ErrL2 = sqrt(semil2) / sqrt(NormExL2);
ErrH1 = sqrt(semih1) / sqrt(NormExH1);

return