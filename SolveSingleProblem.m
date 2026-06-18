clear; addpath src;

%Problem parameters
param.K=6; %wavenumber
param.L=2*pi; %period
param.theta=-pi/5; %incident angle
param.alp=param.K*cos(param.theta); %quasi-periodicity parameter

%Discretizazion parameters
param.h=0.75; %mesh width
param.nd=20; %number of directions
param.M = GetOrder(param); %number of Fourier modes
param.alpha=1/2; param.beta=1/2; param.delta=1/2; %TDG flux coefficients

%define the plane wave direction vectors
param.d=zeros(2,param.nd); 
for l=1:param.nd
    param.d(:,l)=[cos((2*pi*l)/param.nd); sin((2*pi*l)/param.nd)];
end

%% Define domain properties and generate mesh for all layers
NLayers=4; %number of layers

domain1.Hp=12; domain1.Hm=8; domain1.epsilon=[1, 2]; domain1.name='double_rectangle';
[Tmat{1}.mesh,Tmat{1}.param] = GenerateMesh(param,domain1);
Tmat{1}.Hp=domain1.Hp; Tmat{1}.Hm=domain1.Hm; 

domain2.Hp=6; domain2.Hm=2; domain2.epsilon=2; domain2.name='dir_poly'; 
domain2.vertices=[3, 0.5; 5, 0.5; 5, -0.5; 3, -0.5];
[Tmat{2}.mesh,Tmat{2}.param] = GenerateMesh(param,domain2); 
Tmat{2}.Hp=domain2.Hp; Tmat{2}.Hm=domain2.Hm; 

domain3.Hp=0; domain3.Hm=-4; domain3.epsilon=[2, 1.5]; domain3.name='u_shape';
[Tmat{3}.mesh,Tmat{3}.param] = GenerateMesh(param,domain3);
Tmat{3}.Hp=domain3.Hp; Tmat{3}.Hm=domain3.Hm; 

domain4.Hp=-8; domain4.Hm=-12; domain4.epsilon=[1.5, (1.27+0.1*1i)^2, 1.5]; 
domain4.name='penetrable_poly'; domain4.vertices=[2, 0; 3, 1; 4, 0; 3, -1];
[Tmat{4}.mesh,Tmat{4}.param] = GenerateMesh(param,domain4);
Tmat{4}.Hp=domain4.Hp; Tmat{4}.Hm=domain4.Hm; 

%% Compute T-matrices in parallel
parfor j = 1:NLayers
    Tmat{j} = BuildTMatrix(Tmat{j});
end

%% Build and solve coupled linear system
M=param.M; Ndim=2*M+1;

%incident vector
a_d=zeros(Ndim,1); a_d(M+1)=exp(-1i*Tmat{1}.Hp*Tmat{1}.beta_n_p(M+1));
a_u=zeros(Ndim,1);

s = MultiLayerSolve(NLayers,Tmat,a_u,a_d);

%% Plot solution between layers
C=2; %parameter for limit above and below layers
MultiLayerPlot(Tmat,NLayers,s,a_d,a_u,C);

%% Plot inside the layers (if needed)
PlotInsideLayers(Tmat,NLayers,s,a_d,a_u)

%% Plot layers profile and limits
PlotLayerProfile(Tmat,C)



