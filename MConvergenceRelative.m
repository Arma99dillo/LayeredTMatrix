clear; addpath quadtriangle; addpath src

%Problem parameters
param.L=2*pi; %period
param.theta=-pi/3; %incident angle
NLayers=2; %number of layers

%Discretization parameters
param.h=0.75; %mesh width
param.nd=20; %number of directions
M_raff=100; %number of Fourier modes for the refined solution
param.alpha=1/2; param.beta=1/2; param.delta=1/2; %TDG flux coefficients 

%define the plane wave direction vectors
param.d=zeros(2,param.nd);
for l=1:param.nd
    param.d(:,l)=[cos((2*pi*l)/param.nd); sin((2*pi*l)/param.nd)];
end

%% Define domain properties and generate mesh for both layers
domain1.Hp=5; domain1.Hm=1; domain1.epsilon=[1, 3]; domain1.name='u_shape';
[Tmat{1}.mesh,Tmat{1}.param] = GenerateMesh(param,domain1);
Tmat{1}.Hp=domain1.Hp; Tmat{1}.Hm=domain1.Hm; 

domain2.Hp=0; domain2.Hm=-4; domain2.epsilon=[3, (1.27+0.1*1i)^2, 3]; 
domain2.name='penetrable_poly'; domain2.vertices=[2, 0; 3, 1; 4, 0; 3, -1];
[Tmat{2}.mesh,Tmat{2}.param] = GenerateMesh(param,domain2);
Tmat{2}.Hp=domain2.Hp; Tmat{2}.Hm=domain2.Hm; 

%% M-convergence test
KK=[4,6,8,10]; %different wavenumbers
Mi=2; Mf=50; %min and max value of M
L2Error=zeros((Mf-Mi)/2+1,size(KK,2)); %error vectors

%cycle on K and M
for j=1:size(KK,2)

    %get wavenumber and quasi-periodicity parameter
    Tmat{1}.param.K=KK(j); Tmat{2}.param.K=KK(j);
    Tmat{1}.param.alp=Tmat{1}.param.K*cos(param.theta); 
    Tmat{2}.param.alp=Tmat{1}.param.alp; 

    disp(['M-convergence test, k=', num2str(KK(j))])

    %Refined solution for relative error
    Tmat_raff = Tmat;
    Tmat_raff{1}.param.M=M_raff;
    Tmat_raff{2}.param.M=M_raff;

    disp(['Computing refined solution (M=', num2str(M_raff), ')' ])

    %Build T-matrices for the refined solution
    Tmat_raff{1} = BuildTMatrix(Tmat_raff{1});
    Tmat_raff{2} = BuildTMatrix(Tmat_raff{2});

    %Build and solve coupled linear system
    Ndim_raff=2*M_raff+1; nsol=0;
    a_d_raff=zeros(Ndim_raff,1); a_d_raff(nsol+M_raff+1)=exp(-1i*Tmat_raff{1}.Hp*Tmat_raff{1}.beta_n_p(nsol+M_raff+1));
    a_u_raff=zeros(Ndim_raff,1);

    s_raff = MultiLayerSolve(NLayers,Tmat_raff,a_u_raff,a_d_raff);

    disp(['Computed refined solution (M=', num2str(M_raff), ')' ])

    %Cycle on number of Fourier modes
    v=1;
    for M=Mi:2:Mf

        Tmat{1}.param.M=M; Tmat{2}.param.M=M;

        %Build T-matrices 
        Tmat{1} = BuildTMatrix(Tmat{1}); Tmat{2} = BuildTMatrix(Tmat{2});

        %Build and solve coupled linear system
        Ndim=2*M+1; nsol=0;
        a_d=zeros(Ndim,1); a_d(nsol+M+1)=exp(-1i*Tmat{1}.Hp*Tmat{1}.beta_n_p(nsol+M+1));
        a_u=zeros(Ndim,1);
    
        s = MultiLayerSolve(NLayers,Tmat,a_u,a_d);

        %Compute error against refined solution
        [err2,~] = MultiLayerErrRel(Tmat,Tmat_raff,NLayers,s,s_raff,a_d,a_u,a_d_raff,a_u_raff);
        L2Error(v,j) = err2;
        disp([ 'Computed error for k=', num2str(KK(j)) , ', M=', num2str(M), ' Fourier modes truncation' ])

        v=v+1;
    end

end

%% Convergence plot
figure()
ax=gca;
semilogy(Mi:2:Mf,L2Error(:,1),'o-',Mi:2:Mf,L2Error(:,2),'*-', ...
    Mi:2:Mf,L2Error(:,3),'^-',Mi:2:Mf,L2Error(:,4),'s-','LineWidth',1.2);
grid on; ax.GridAlpha = 0.15;
xline(0.5.*KK.*param.L,'--k',{'k=4','k=6','k=8','k=10'},'LineWidth',1.2,'FontSize',14);
xlim([Mi Mf]); xticks(5:5:Mf); 
LL = legend('$k=4$','$k=6$','$k=8$','$k=10$','FontSize', 14);
set(LL, 'Interpreter', 'latex');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
xlabel('Number of Fourier modes','FontSize',18, 'Interpreter','latex')
ylabel('$L^2$ Error','FontSize',18, 'Interpreter','latex')

%% Save plot
exportgraphics(gcf, 'figs/MConv.pdf', 'ContentType', 'vector')