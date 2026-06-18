function MultiLayerPlotSameType(LayerArr,Tmat,b,a_d,a_u,C)

%plots approximate solution ouside the layers in [0,L]x[Hm-C,Hm+C]

%get parameters
M=Tmat{1}.param.M; L=Tmat{1}.param.L; alpha_n=Tmat{1}.alpha_n;
Ndim=2*M+1;
shape=LayerArr.shape; pos=LayerArr.pos; refl=LayerArr.refl;
NLayer=size(shape,1);

figure()
disp('Computing and plotting the T-matrix solution between the layers')
tic()

%functions Phi_n to evaluate the T-matrix solution
eval_up = @(x1,x2,alpha,beta,Hp) exp(1i.*beta.*(x2-Hp)).*exp(1i.*alpha.*x1);
eval_down = @(x1,x2,alpha,beta,Hm) exp(-1i.*beta.*(x2-Hm)).*exp(1i.*alpha.*x1);

%setup grid of vertical width C on the upper side
npoints=500;
g1=linspace(0,L,npoints); g2=linspace(Tmat{shape(1)}.Hp+pos(1),Tmat{shape(1)}.Hp+pos(1)+C,npoints);
[x1,x2]=meshgrid(g1,g2);

%get vector for upward and downward expansion and compute solution 
%using the T-matrix expansion
b_u=b(1:Ndim);

%reflect T-matrix if needed
if refl(1) == 1
    Tmat_temp=ReflectTmat(Tmat{shape(1)});
else
    Tmat_temp=Tmat{shape(1)};
end 

sol=zeros(npoints,npoints); v=1;
for n=-M:M
    sol=sol+b_u(v).*eval_up(x1,x2,alpha_n(v),Tmat_temp.beta_n_p(v),Tmat_temp.Hp+pos(1));
    sol=sol+a_d(v).*eval_down(x1,x2,alpha_n(v),Tmat_temp.beta_n_p(v),Tmat_temp.Hp+pos(1));
    v=v+1;
end

%plot solution
subplot(1,2,1)
surf(x1,x2,real(sol)); grid off; shading interp; colorbar
view(2); 
hold on

subplot(1,2,2)
surf(x1,x2,abs(sol)); grid off; shading interp; colorbar
view(2); 
hold on

r=1;
for j=1:NLayer-1

    %setup grid between layers
    npoints=200;
    g1=linspace(0,L,npoints); g2=linspace(Tmat{shape(j+1)}.Hp+pos(j+1),Tmat{shape(j)}.Hm+pos(j),npoints);
    [x1,x2]=meshgrid(g1,g2);
    
    %get vector for upward and downward expansion and compute solution 
%using the T-matrix expansion
    b_d=b(r*Ndim+1:(r+1)*Ndim); b_u=b((r+1)*Ndim+1:(r+2)*Ndim); 

    if refl(j) == 1
        Tmat_temp1=ReflectTmat(Tmat{shape(j)});
    else
        Tmat_temp1=Tmat{shape(j)};
    end

    if refl(j+1) == 1
        Tmat_temp2=ReflectTmat(Tmat{shape(j)});
    else
        Tmat_temp2=Tmat{shape(j)};
    end

    sol=zeros(npoints,npoints); v=1;
    for n=-M:M
        sol=sol+b_u(v).*eval_up(x1,x2,alpha_n(v),Tmat_temp2.beta_n_p(v),Tmat_temp2.Hp+pos(j+1));
        sol=sol+b_d(v).*eval_down(x1,x2,alpha_n(v),Tmat_temp1.beta_n_m(v),Tmat_temp1.Hm+pos(j));
        v=v+1;
    end

    %plot solution
    subplot(1,2,1)
    surf(x1,x2,real(sol)); grid off; shading interp; colorbar
    view(2);
    hold on

    subplot(1,2,2)
    surf(x1,x2,abs(sol)); grid off; shading interp; colorbar
    view(2);
    hold on

    r=r+2;
end

%setup a grid of vertical width C on the lower side
g1=linspace(0,L,npoints); g2=linspace(Tmat{shape(NLayer)}.Hm+pos(NLayer)-C,Tmat{shape(NLayer)}.Hm+pos(NLayer),npoints);
[x1,x2]=meshgrid(g1,g2);

%get vector for upward and downward expansion and compute solution 
%using the T-matrix expansion
b_d=b(2*NLayer*Ndim-Ndim+1:end);

if refl(NLayer) == 1
    Tmat_temp=ReflectTmat(Tmat{shape(NLayer)});
else
    Tmat_temp=Tmat{shape(NLayer)};
end

sol=zeros(npoints,npoints); v=1;
for n=-M:M
    sol=sol+a_u(v).*eval_up(x1,x2,alpha_n(v),Tmat_temp.beta_n_m(v),Tmat_temp.Hm+pos(NLayer));
    sol=sol+b_d(v).*eval_down(x1,x2,alpha_n(v),Tmat_temp.beta_n_m(v),Tmat_temp.Hm+pos(NLayer));
    v=v+1;
end

%plot solution
subplot(1,2,1)
surf(x1,x2,real(sol)); grid off; shading interp; colorbar
view(2); 
hold on

subplot(1,2,2)
surf(x1,x2,abs(sol)); grid off; shading interp; colorbar
view(2); 
hold on

toc()

return