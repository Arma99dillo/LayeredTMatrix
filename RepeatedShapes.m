clear; addpath src;

%Problem parameters
param.K=8; %wavenumber
param.L=4; %period
param.theta=-pi/4; %incident angle
param.alp=param.K*cos(param.theta); %quasi-periodicity parameter

%Discretizazion parameters
param.h=0.5; %mesh width
param.nd=20; %number of directions
param.M=GetOrder(param); %number of Fourier modes
param.alpha=1/2; param.beta=1/2; param.delta=1/2; %TDG flux coefficients

%define the plane wave direction vectors
param.d=zeros(2,param.nd); 
for l=1:param.nd
    param.d(:,l)=[cos((2*pi*l)/param.nd); sin((2*pi*l)/param.nd)];
end

%% Define domain properties and generate mesh for all layers shape
NShape=2; %number of shapes
LayerShape=cell(NShape,1);

%First shape
LayerShape{1}.Hp=1; LayerShape{1}.Hm=-1; LayerShape{1}.epsilon=[1, 2, 1]; 
LayerShape{1}.name='penetrable_poly'; LayerShape{1}.vertices=[2,-0.5;2.5,0.5;3,-0.5];

%Second shape
LayerShape{2}.Hp=1; LayerShape{2}.Hm=-1; LayerShape{2}.epsilon=1; 
LayerShape{2}.name='dir_poly'; LayerShape{2}.vertices=[1, 0; 1.5, 0.5; 2, 0; 1.5, -0.5];

%Generate mesh for every shape in parallel
parfor j=1:NShape
    [Tmat{j}.mesh,Tmat{j}.param] = GenerateMesh(param,LayerShape{j});
    Tmat{j}.Hp=LayerShape{j}.Hp; Tmat{j}.Hm=LayerShape{j}.Hm;
    Tmat{j}.shape=LayerShape{j};
end

%% Compute T-matrix for every shape in parallel
parfor j=1:NShape
    disp(['T-Matrix ', num2str(j), ' computation'])
    Tmat{j} = BuildTMatrix(Tmat{j});
end

%% Build coupled linear system
%structure containing the layers arrangement properties
LayerArr.shape = [1; 2; 1; 1; 2; 1; 2; 2; 1; 2]; %shape
LayerArr.pos = [11;9;7;3;1;-2;-4;-6;-9;-11]; %position of the center
LayerArr.refl = [0;0;1;0;0;1;0;0;0;0]; %reflection

%incident vector
M=Tmat{1}.param.M; Ndim=2*M+1;
a_d=zeros(Ndim,1); a_d(M+1)=exp(-1i*(Tmat{1}.Hp+LayerArr.pos(1))*Tmat{1}.beta_n_p(M+1));
a_u=zeros(Ndim,1);

%build and solve coupled system
s = MultiLayerSolveSameType(LayerArr,Tmat,a_u,a_d);

%% Plot solution between layers
C=2; %parameter for limit above and below layers
MultiLayerPlotSameType(LayerArr,Tmat,s,a_d,a_u,C);

%% Plot inside the layers (if needed)
PlotInsideLayersSameType(LayerArr,Tmat,s,a_d,a_u)

%% Plot layers profile and limits
PlotLayerProfileSameType(LayerArr,Tmat,C);