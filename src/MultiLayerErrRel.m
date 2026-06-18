function [ErrL2,ErrH1] = MultiLayerErrRel(Tmat,Tmat_raff,Nlayer,b,b_raff,a_d,a_u,a_d_raff,a_u_raff)

%compute relative error against refined solution

%get parameters
M=Tmat{1}.param.M; L=Tmat{1}.param.L; alpha_n=Tmat{1}.alpha_n; Ndim=2*M+1;
M_raff=Tmat_raff{1}.param.M; Ndim_raff=2*M_raff+1; alpha_n_raff=Tmat_raff{1}.alpha_n;

disp('Computing error against refined solution')
tic()

%functions Phi_n to evaluate the T-matrix solution
eval_up = @(x1,x2,alpha,beta,Hp) exp(1i.*beta.*(x2-Hp)).*exp(1i.*alpha.*x1);
eval_down = @(x1,x2,alpha,beta,Hm) exp(-1i.*beta.*(x2-Hm)).*exp(1i.*alpha.*x1);

semil2=0; semih1=0; NormExL2=0; NormExH1=0;

%first compute error between the layers
r=1;
for j=1:Nlayer-1

    %setup grid in the middle
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
    
    %compute refined solution using the T-matrix expansion
    b_d_r=b_raff(r*Ndim_raff+1:(r+1)*Ndim_raff); b_u_r=b_raff((r+1)*Ndim_raff+1:(r+2)*Ndim_raff); 
    sol_raff=zeros(npoints,npoints); g1_raff=zeros(npoints,npoints); g2_raff=zeros(npoints,npoints); v=1;
    for n=-M_raff:M_raff
        sol_raff=sol_raff+b_u_r(v).*eval_up(x1,x2,alpha_n_raff(v),Tmat_raff{j+1}.beta_n_p(v),Tmat_raff{j+1}.Hp);
        sol_raff=sol_raff+b_d_r(v).*eval_down(x1,x2,alpha_n_raff(v),Tmat_raff{j}.beta_n_m(v),Tmat_raff{j}.Hm);
        g1_raff=g1_raff+1i.*alpha_n_raff(v).*b_u_r(v).*eval_up(x1,x2,alpha_n_raff(v),Tmat_raff{j+1}.beta_n_p(v),Tmat_raff{j+1}.Hp);
        g1_raff=g1_raff+1i.*alpha_n_raff(v).*b_d_r(v).*eval_down(x1,x2,alpha_n_raff(v),Tmat_raff{j}.beta_n_m(v),Tmat_raff{j}.Hm);
        g2_raff=g2_raff+1i.*Tmat_raff{j+1}.beta_n_p(v).*b_u_r(v).*eval_up(x1,x2,alpha_n_raff(v),Tmat_raff{j+1}.beta_n_p(v),Tmat_raff{j+1}.Hp);
        g2_raff=g2_raff-1i.*Tmat_raff{j}.beta_n_m(v).*b_d_r(v).*eval_down(x1,x2,alpha_n_raff(v),Tmat_raff{j}.beta_n_m(v),Tmat_raff{j}.Hm);
        v=v+1;
    end

    %compute error against refined solution
    semil2 = semil2 + sum(abs(sol - sol_raff).^2, 'all') * dx * dy;
    semih1 = semih1 + sum(abs(sol - sol_raff).^2, 'all') * dx * dy + sum(abs(g1 - g1_raff).^2, 'all') * dx * dy + sum(abs(g2 - g2_raff).^2, 'all') * dx * dy;
    NormExL2 = NormExL2 + sum(abs(sol_raff).^2, 'all') * dx * dy;
    NormExH1 = NormExH1 + sum(abs(sol_raff).^2, 'all') * dx * dy + sum(abs(g1_raff).^2, 'all') * dx * dy + sum(abs(g2_raff).^2, 'all') * dx * dy;

    r=r+2;
end

%then compute the error inside layers

%first layer: translate mesh
mesh_temp=Tmat{1}.mesh;
mesh_temp.p(:,2) = mesh_temp.p(:,2)+(Tmat{1}.Hp-Tmat{1}.param.H); %translate mesh
deltaH_p=0; deltaH_m=Tmat{1}.Hm-Tmat{2}.Hp;

%get T-matrix expansion vectors
b_u_p=b(1:Ndim); b_d_p=a_d; b_d_m=b(Ndim+1:2*Ndim); b_u_m=b(2*Ndim+1:3*Ndim);
b_u_p_raff=b_raff(1:Ndim_raff); b_d_p_raff=a_d_raff; b_d_m_raff=b_raff(Ndim_raff+1:2*Ndim_raff); b_u_m_raff=b_raff(2*Ndim_raff+1:3*Ndim_raff);

%solve impedance boundary value problem for both solutions
A = MatrixTDGImp(mesh_temp,Tmat{1}.param);
rhs = rhsTDGImp(mesh_temp,Tmat{1},b_u_p,b_d_p,b_u_m,b_d_m,deltaH_p,deltaH_m);
A_raff = MatrixTDGImp(mesh_temp,Tmat_raff{1}.param);
rhs_raff = rhsTDGImp(mesh_temp,Tmat_raff{1},b_u_p_raff,b_d_p_raff,b_u_m_raff,b_d_m_raff,deltaH_p,deltaH_m);
u = A\rhs; u_raff=A_raff\rhs_raff;

%compute error against refined solution
[errl2,norml2,errh1,normh1]=InsideErrorRel(mesh_temp,Tmat{1}.param,u,Tmat_raff{1}.param,u_raff);
semil2=semil2+errl2;
NormExL2=NormExL2+norml2;
semih1=semih1+errh1;
NormExH1=NormExH1+normh1;


%cycle on middle layers
r=1;
for j=2:Nlayer-1
    mesh_temp=Tmat{j}.mesh;
    mesh_temp.p(:,2) = mesh_temp.p(:,2)+(Tmat{j}.Hp-Tmat{j}.param.H); %translate mesh
    deltaH_p=Tmat{j-1}.Hm-Tmat{j}.Hp; deltaH_m=Tmat{j}.Hm-Tmat{j+1}.Hp;

    %get T-matrix expansion vectors
    b_d_p=b(r*Ndim+1:(r+1)*Ndim); b_u_p=b((r+1)*Ndim+1:(r+2)*Ndim); 
    b_d_m=b((r+2)*Ndim+1:(r+3)*Ndim); b_u_m=b((r+3)*Ndim+1:(r+4)*Ndim);
    b_d_p_raff=b_raff(r*Ndim_raff+1:(r+1)*Ndim_raff); b_u_p_raff=b_raff((r+1)*Ndim_raff+1:(r+2)*Ndim_raff); 
    b_d_m_raff=b_raff((r+2)*Ndim_raff+1:(r+3)*Ndim_raff); b_u_m_raff=b_raff((r+3)*Ndim_raff+1:(r+4)*Ndim_raff);

    %solve impedance boundary value problem for both solutions
    A = MatrixTDGImp(mesh_temp,Tmat{j}.param);
    rhs = rhsTDGImp(mesh_temp,Tmat{j},b_u_p,b_d_p,b_u_m,b_d_m,deltaH_p,deltaH_m);
    A_raff = MatrixTDGImp(mesh_temp,Tmat_raff{j}.param);
    rhs_raff = rhsTDGImp(mesh_temp,Tmat_raff{j},b_u_p_raff,b_d_p_raff,b_u_m_raff,b_d_m_raff,deltaH_p,deltaH_m);
    u = A\rhs; u_raff=A_raff\rhs_raff;

    %compute error against refined solution
    [errl2,norml2,errh1,normh1]=InsideErrorRel(mesh_temp,Tmat{j}.param,u,Tmat_raff{j}.param,u_raff);
    semil2=semil2+errl2;
    NormExL2=NormExL2+norml2;
    semih1=semih1+errh1;
    NormExH1=NormExH1+normh1;

    r=r+2;
end


%last layer: translate mesh
mesh_temp=Tmat{Nlayer}.mesh;
mesh_temp.p(:,2) = mesh_temp.p(:,2)+(Tmat{Nlayer}.Hp-Tmat{Nlayer}.param.H); %translate mesh
deltaH_p=Tmat{Nlayer-1}.Hm-Tmat{Nlayer}.Hp; deltaH_m=0;

%get T-matrix expansion vectors
b_d_p=b(2*Nlayer*Ndim-3*Ndim+1:2*Nlayer*Ndim-2*Ndim); b_u_p=b(2*Nlayer*Ndim-2*Ndim+1:2*Nlayer*Ndim-Ndim);
b_d_m=b(2*Nlayer*Ndim-Ndim+1:end); b_u_m=a_u;
b_d_p_raff=b_raff(2*Nlayer*Ndim_raff-3*Ndim_raff+1:2*Nlayer*Ndim_raff-2*Ndim_raff); b_u_p_raff=b_raff(2*Nlayer*Ndim_raff-2*Ndim_raff+1:2*Nlayer*Ndim_raff-Ndim_raff);
b_d_m_raff=b_raff(2*Nlayer*Ndim_raff-Ndim_raff+1:end); b_u_m_raff=a_u_raff;

%solve impedance boundary value problem for both solutions
A = MatrixTDGImp(mesh_temp,Tmat{Nlayer}.param);
rhs = rhsTDGImp(mesh_temp,Tmat{Nlayer},b_u_p,b_d_p,b_u_m,b_d_m,deltaH_p,deltaH_m);
A_raff = MatrixTDGImp(mesh_temp,Tmat_raff{Nlayer}.param);
rhs_raff = rhsTDGImp(mesh_temp,Tmat_raff{Nlayer},b_u_p_raff,b_d_p_raff,b_u_m_raff,b_d_m_raff,deltaH_p,deltaH_m);
u = A\rhs; u_raff=A_raff\rhs_raff;

%compute error against refined solution
[errl2,norml2,errh1,normh1]=InsideErrorRel(mesh_temp,Tmat{Nlayer}.param,u,Tmat_raff{Nlayer}.param,u_raff);
semil2=semil2+errl2;
NormExL2=NormExL2+norml2;
semih1=semih1+errh1;
NormExH1=NormExH1+normh1;

toc()

%return relative errors
ErrL2 = sqrt(semil2) / sqrt(NormExL2);
ErrH1 = sqrt(semih1) / sqrt(NormExH1);

return