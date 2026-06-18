function PlotInsideLayersSameType(LayerArr,Tmat,b,a_d,a_u)

%plots approximate solution inside the layers

%get parameters
M=Tmat{1}.param.M; Ndim=2*M+1;
shape=LayerArr.shape; pos=LayerArr.pos; refl=LayerArr.refl;
NLayer=size(shape,1);

%first layer: translate mesh
disp('Plot inside layer 1')
mesh_temp=Tmat{shape(1)}.mesh;
mesh_temp.p(:,2) = mesh_temp.p(:,2)+pos(1); %translate mesh

%get T-matrix expansion vectors
b_u_p=b(1:Ndim); b_d_p=a_d; b_d_m=b(Ndim+1:2*Ndim); b_u_m=b(2*Ndim+1:3*Ndim);
deltaH_p=0; deltaH_m=Tmat{shape(1)}.Hm+pos(1)-Tmat{shape(2)}.Hp-pos(2);
%reflect T-matrix if needed
if refl(1) == 1
    Tmat_temp=ReflectTmat(Tmat{shape(1)});
else
    Tmat_temp=Tmat{shape(1)};
end 
Tmat_temp.Hp=Tmat{shape(1)}.Hp+pos(1); Tmat_temp.Hm=Tmat{shape(1)}.Hm+pos(1);

%solve impedance boundary value problem
A = MatrixTDGImp(mesh_temp,Tmat_temp.param);
rhs = rhsTDGImp(mesh_temp,Tmat_temp,b_u_p,b_d_p,b_u_m,b_d_m,deltaH_p,deltaH_m);
u = A\rhs;

%plot solution
PlotSolution(mesh_temp,Tmat_temp.param,u);
hold on

%cycle on middle layers
r=1;
for j=2:NLayer-1

    %translate mesh
    disp(['Plot inside layer ', num2str(j)])
    mesh_temp=Tmat{shape(j)}.mesh;
    mesh_temp.p(:,2) = mesh_temp.p(:,2)+pos(j); %translate mesh

    %get T-matrix expansion vectors
    b_d_p=b(r*Ndim+1:(r+1)*Ndim); b_u_p=b((r+1)*Ndim+1:(r+2)*Ndim); 
    b_d_m=b((r+2)*Ndim+1:(r+3)*Ndim); b_u_m=b((r+3)*Ndim+1:(r+4)*Ndim);
    deltaH_p=Tmat{shape(j-1)}.Hm+pos(j-1)-Tmat{shape(j)}.Hp-pos(j);
    deltaH_m=Tmat{shape(j)}.Hm+pos(j)-Tmat{shape(j+1)}.Hp-pos(j+1);
    if refl(j) == 1
        Tmat_temp=ReflectTmat(Tmat{shape(j)});
    else
        Tmat_temp=Tmat{shape(j)};
    end
    Tmat_temp.Hp=Tmat{shape(j)}.Hp+pos(j); Tmat_temp.Hm=Tmat{shape(j)}.Hm+pos(j);

    %solve impedance boundary value problem
    A = MatrixTDGImp(mesh_temp,Tmat_temp.param);
    rhs = rhsTDGImp(mesh_temp,Tmat_temp,b_u_p,b_d_p,b_u_m,b_d_m,deltaH_p,deltaH_m);
    u = A\rhs;

    %plot solution
    PlotSolution(mesh_temp,Tmat_temp.param,u);
    hold on

    r=r+2;
end


%last layer:translate mesh
disp(['Plot inside layer ', num2str(NLayer)])
mesh_temp=Tmat{shape(NLayer)}.mesh;
mesh_temp.p(:,2) = mesh_temp.p(:,2)+pos(NLayer); %translate mesh

%get T-matrix expansion vectors
b_d_p=b(2*NLayer*Ndim-3*Ndim+1:2*NLayer*Ndim-2*Ndim); b_u_p=b(2*NLayer*Ndim-2*Ndim+1:2*NLayer*Ndim-Ndim);
b_d_m=b(2*NLayer*Ndim-Ndim+1:end); b_u_m=a_u;
deltaH_p=Tmat{shape(NLayer-1)}.Hm+pos(NLayer-1)-Tmat{shape(NLayer)}.Hp-pos(NLayer); deltaH_m=0;
if refl(NLayer) == 1
    Tmat_temp=ReflectTmat(Tmat{shape(NLayer)});
else
    Tmat_temp=Tmat{shape(NLayer)};
end
Tmat_temp.Hp=Tmat{shape(NLayer)}.Hp+pos(NLayer); Tmat_temp.Hm=Tmat{shape(NLayer)}.Hm+pos(NLayer);

%solve impedance boundary value problem
A = MatrixTDGImp(mesh_temp,Tmat_temp.param);
rhs = rhsTDGImp(mesh_temp,Tmat_temp,b_u_p,b_d_p,b_u_m,b_d_m,deltaH_p,deltaH_m);
u = A\rhs;

%plot solution
PlotSolution(mesh_temp,Tmat_temp.param,u);
hold on

end