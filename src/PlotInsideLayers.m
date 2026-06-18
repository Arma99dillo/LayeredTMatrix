function PlotInsideLayers(Tmat,Nlayer,b,a_d,a_u)

%plots approximate solution inside the layers

%get parameters
M=Tmat{1}.param.M; Ndim=2*M+1;

%first layer: translate mesh
disp('Plot inside layer 1')
mesh_temp=Tmat{1}.mesh;
mesh_temp.p(:,2) = mesh_temp.p(:,2)+(Tmat{1}.Hp-Tmat{1}.param.H); %translate mesh

%get T-matrix expansion vectors
b_u_p=b(1:Ndim); b_d_p=a_d; b_d_m=b(Ndim+1:2*Ndim); b_u_m=b(2*Ndim+1:3*Ndim);
deltaH_p=0; deltaH_m=Tmat{1}.Hm-Tmat{2}.Hp;

%solve impedance boundary value problem
A = MatrixTDGImp(mesh_temp,Tmat{1}.param);
rhs = rhsTDGImp(mesh_temp,Tmat{1},b_u_p,b_d_p,b_u_m,b_d_m,deltaH_p,deltaH_m);
u = A\rhs;

%plot solution
PlotSolution(mesh_temp,Tmat{1}.param,u);
hold on

%cycle on middle layers
r=1;
for j=2:Nlayer-1

    %translate mesh
    disp(['Plot inside layer ', num2str(j)])
    mesh_temp=Tmat{j}.mesh;
    mesh_temp.p(:,2) = mesh_temp.p(:,2)+(Tmat{j}.Hp-Tmat{j}.param.H); %translate mesh

    %get T-matrix expansion vectors
    b_d_p=b(r*Ndim+1:(r+1)*Ndim); b_u_p=b((r+1)*Ndim+1:(r+2)*Ndim); 
    b_d_m=b((r+2)*Ndim+1:(r+3)*Ndim); b_u_m=b((r+3)*Ndim+1:(r+4)*Ndim);
    deltaH_p=Tmat{j-1}.Hm-Tmat{j}.Hp; deltaH_m=Tmat{j}.Hm-Tmat{j+1}.Hp;

    %solve impedance boundary value problem
    A = MatrixTDGImp(mesh_temp,Tmat{j}.param);
    rhs = rhsTDGImp(mesh_temp,Tmat{j},b_u_p,b_d_p,b_u_m,b_d_m,deltaH_p,deltaH_m);
    u = A\rhs;

    %plot solution
    PlotSolution(mesh_temp,Tmat{j}.param,u);
    hold on
    r=r+2;
end

%last layer: translate mesh
disp(['Plot inside layer ', num2str(Nlayer)])
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

%plot solution
PlotSolution(mesh_temp,Tmat{Nlayer}.param,u);
hold on


end