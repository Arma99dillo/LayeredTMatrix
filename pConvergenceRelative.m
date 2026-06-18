clear; addpath src; addpath quadtriangle;

%Problem parameters
param.K=4; %wavenumber
param.L=2*pi; %period
param.theta=-pi/3; %incident angle
param.alp=param.K*cos(param.theta); %quasi-periodicity parameter
NLayers=2; %number of layers

%Discretizazion parameters
param.h=0.75; %mesh width
param.M = GetOrder(param); %number of Fourier modes
param.alpha=1/2; param.beta=1/2; param.delta=1/2; %TDG flux coefficients 

%% Define domain properties and generate mesh for both layers
domain1.Hp=5; domain1.Hm=1; %Height delimiting the layer
domain1.epsilon=[1, 3]; domain1.name='u_shape'; %value of epsilon and shape
[Tmat{1}.mesh,Tmat{1}.param] = GenerateMesh(param,domain1); %generate mesh
Tmat{1}.Hp=domain1.Hp; Tmat{1}.Hm=domain1.Hm; %assing value of Hp and Hm

domain2.Hp=0; domain2.Hm=-4;  %Height delimiting the layer
domain2.epsilon=[3, (1.27+0.1*1i)^2, 3]; %value of epsilon and shape
domain2.name='penetrable_poly'; domain2.vertices=[2, 0; 3, 1; 4, 0; 3, -1];
[Tmat{2}.mesh,Tmat{2}.param] = GenerateMesh(param,domain2); %generate mesh
Tmat{2}.Hp=domain2.Hp; Tmat{2}.Hm=domain2.Hm; %assing value of Hp and Hm

%% Compute refined solution
p_ref=20; %number of directions for the refined solution

disp('Computing refined solution')

Tmat_raff = Tmat;
Tmat_raff{1}.param.nd=p_ref; Tmat_raff{2}.param.nd=p_ref;
Tmat_raff{1}.param.d=zeros(2,Tmat_raff{1}.param.nd);
for l=1:Tmat_raff{1}.param.nd %define equispaced directions
    Tmat_raff{1}.param.d(:,l)=[cos((2*pi*l)/Tmat_raff{1}.param.nd); sin((2*pi*l)/Tmat_raff{1}.param.nd)];
end
Tmat_raff{2}.param.d=Tmat_raff{1}.param.d;

%compute T-matrices
Tmat_raff{1} = BuildTMatrix(Tmat_raff{1}); 
Tmat_raff{2} = BuildTMatrix(Tmat_raff{2});

%Build and solve coupling system
M_raff=Tmat_raff{1}.param.M; Ndim_raff=2*M_raff+1; 
a_d_raff=zeros(Ndim_raff,1); 
a_d_raff(M_raff+1)=exp(-1i*Tmat_raff{1}.Hp*Tmat_raff{1}.beta_n_p(M_raff+1)); %incident vector from above
a_u_raff=zeros(Ndim_raff,1);

s_raff = MultiLayerSolve(NLayers,Tmat_raff,a_u_raff,a_d_raff); 


%% Cycle on number of directions
ni=3; nf=p_ref-1; %min and max value of directions
L2Error=zeros(nf-ni+1,1); H1Error=zeros(nf-ni+1,1); v=1; %error vector

for p=ni:nf %cycle on number of PW directions

    disp(['p= ',num2str(p)])

    Tmat{1}.param.nd=p; Tmat{2}.param.nd=p; %number of directions
    Tmat{1}.param.d=zeros(2,Tmat{1}.param.nd); 
    for l=1:Tmat{1}.param.nd %define equispaced directions
        Tmat{1}.param.d(:,l)=[cos((2*pi*l)/Tmat{1}.param.nd); sin((2*pi*l)/Tmat{1}.param.nd)];
    end  
    Tmat{2}.param.d=Tmat{1}.param.d;

    %compute T-matrices
    Tmat{1} = BuildTMatrix(Tmat{1}); Tmat{2} = BuildTMatrix(Tmat{2});

    %Build and solve coupling system
    M=Tmat{1}.param.M; Ndim=2*M+1; 
    a_d=zeros(Ndim,1); a_d(M+1)=exp(-1i*Tmat{1}.Hp*Tmat{1}.beta_n_p(M+1)); %incident vector from above
    a_u=zeros(Ndim,1);

    s = MultiLayerSolve(NLayers,Tmat,a_u,a_d);

    %Compute relative error against refined solution
    [err2,err1] = MultiLayerErrRel(Tmat,Tmat_raff,NLayers,s,s_raff,a_d,a_u,a_d_raff,a_u_raff);
    L2Error(v) = err2; H1Error(v) = err1;
    v=v+1;
end

%% Plot solution and error
MultiLayerPlotErrRel(Tmat,Tmat_raff,NLayers,s,s_raff,a_d,a_u,a_d_raff,a_u_raff)

%% Plot layers profile and limits
C=0; %parameter for limit above and below layers
PlotLayerProfile(Tmat,C)

%% Plot error convergence
PlotConvergence(ni,nf,L2Error,H1Error,ni:2:nf)

%% Save plot
exportgraphics(gcf, 'figs/RelativeConv.pdf', 'ContentType', 'vector')