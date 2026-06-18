function MultiLayerPlot(Tmat,Nlayer,b,a_d,a_u,C)

%plots approximate solution ouside the layers in [0,L]x[Hm-C,Hm+C]

%get parameters
M=Tmat{1}.param.M; L=Tmat{1}.param.L; alpha_n=Tmat{1}.alpha_n;
Ndim=2*M+1;
npoints=200;

figure()
disp('Computing and plotting the T-matrix solution between the layers')
tic()

%functions Phi_n to evaluate the T-matrix solution
eval_up = @(x1,x2,alpha,beta,Hp) exp(1i.*beta.*(x2-Hp)).*exp(1i.*alpha.*x1);
eval_down = @(x1,x2,alpha,beta,Hm) exp(-1i.*beta.*(x2-Hm)).*exp(1i.*alpha.*x1);

%setup grid of vertical width C on the upper side
g1=linspace(0,L,npoints); g2=linspace(Tmat{1}.Hp,Tmat{1}.Hp+C,npoints);
[x1,x2]=meshgrid(g1,g2);

%get vector for upward and downward expansion and compute solution 
%using the T-matrix expansion
b_u=b(1:Ndim);
sol=zeros(npoints,npoints); v=1;
for n=-M:M
    sol=sol+b_u(v).*eval_up(x1,x2,alpha_n(v),Tmat{1}.beta_n_p(v),Tmat{1}.Hp);
    sol=sol+a_d(v).*eval_down(x1,x2,alpha_n(v),Tmat{1}.beta_n_p(v),Tmat{1}.Hp);
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

%cycle on layers
r=1;
for j=1:Nlayer-1

    %setup grid between layers
    g1=linspace(0,L,npoints); g2=linspace(Tmat{j+1}.Hp,Tmat{j}.Hm,npoints);
    [x1,x2]=meshgrid(g1,g2);
    
    %get vector for upward and downward expansion and compute solution 
    %using the T-matrix expansion
    b_d=b(r*Ndim+1:(r+1)*Ndim); b_u=b((r+1)*Ndim+1:(r+2)*Ndim); 
    sol=zeros(npoints,npoints); v=1;
    for n=-M:M
        sol=sol+b_u(v).*eval_up(x1,x2,alpha_n(v),Tmat{j+1}.beta_n_p(v),Tmat{j+1}.Hp);
        sol=sol+b_d(v).*eval_down(x1,x2,alpha_n(v),Tmat{j}.beta_n_m(v),Tmat{j}.Hm);
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

% setup a grid of vertical width C on the lower side
g1=linspace(0,L,npoints); g2=linspace(Tmat{Nlayer}.Hm-C,Tmat{Nlayer}.Hm,npoints);
[x1,x2]=meshgrid(g1,g2);

%get vector for upward and downward expansion and compute solution 
%using the T-matrix expansion
b_d=b(2*Nlayer*Ndim-Ndim+1:end);
sol=zeros(npoints,npoints); v=1;
for n=-M:M
    sol=sol+a_u(v).*eval_up(x1,x2,alpha_n(v),Tmat{Nlayer}.beta_n_m(v),Tmat{Nlayer}.Hm);
    sol=sol+b_d(v).*eval_down(x1,x2,alpha_n(v),Tmat{Nlayer}.beta_n_m(v),Tmat{Nlayer}.Hm);
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