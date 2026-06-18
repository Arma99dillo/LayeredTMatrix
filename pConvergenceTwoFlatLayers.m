clear; addpath src; addpath quadtriangle;

%Problem parameters
param.K=5; %wavenumber
param.L=2*pi; %period
param.theta=-pi/3; %incident angle 
param.alp=param.K*cos(param.theta); %quasi-periodicity parameter
NLayers=2; %number of layers

%Discretization parameters
param.h=1; %mesh width
param.M = GetOrder(param); %number of Fourier modes
param.alpha=1/2; param.beta=1/2; param.delta=1/2; %TDG flux coefficients 

%% Define domain properties and generate mesh for both layers
domain1.Hp=6; domain1.Hm=2; domain1.epsilon=[1, 2]; domain1.name='double_rectangle';
[Tmat{1}.mesh,Tmat{1}.param] = GenerateMesh(param,domain1);
Tmat{1}.Hp=domain1.Hp; Tmat{1}.Hm=domain1.Hm; 

domain2.Hp=-2; domain2.Hm=-6; domain2.epsilon=[2, 1]; domain2.name='double_rectangle';
[Tmat{2}.mesh,Tmat{2}.param] = GenerateMesh(param,domain2);
Tmat{2}.Hp=domain2.Hp; Tmat{2}.Hm=domain2.Hm; 

%% Cycle on number of directions
ni=3; nf=30; %min and max value of directions
L2Error=zeros(nf-ni+1,1); H1Error=zeros(nf-ni+1,1); v=1; %error vectors

for nd=ni:nf %cycle on number of PW directions
    
    disp(['p= ',num2str(nd)])

    Tmat{1}.param.nd=nd; Tmat{2}.param.nd=nd; %number of plane wave directions
    Tmat{1}.param.d=zeros(2,Tmat{1}.param.nd); 
    for l=1:Tmat{1}.param.nd
        Tmat{1}.param.d(:,l)=[cos((2*pi*l)/Tmat{1}.param.nd); sin((2*pi*l)/Tmat{1}.param.nd)];
    end
    Tmat{2}.param.d=Tmat{1}.param.d; 
    
    %compute T-matrices
    Tmat{1} = BuildTMatrix(Tmat{1}); Tmat{2} = BuildTMatrix(Tmat{2});

    %setup incident field expansion vector
    Ndim=2*param.M+1;
    a_d=zeros(Ndim,1); a_d(param.M+1)=exp(-1i*Tmat{1}.Hp*Tmat{1}.beta_n_p(param.M+1));
    a_u=zeros(Ndim,1);

    %Build and solve global linear system
    s = MultiLayerSolve(NLayers,Tmat,a_u,a_d);

    % Error against exact solution
    [err2, err1] = TwoFlatLayersErr(Tmat,NLayers,s,a_d,a_u);
    L2Error(v) = err2; H1Error(v) = err1;
    v=v+1;
end

%% Plot solution and error
TwoFlatLayersPlotErr(Tmat,NLayers,s,a_d,a_u)

%% Plot error
PlotConvergence(ni,nf,L2Error,H1Error,5:5:nf)

%% Save plot
exportgraphics(gcf, 'figs/pConvTwoLayers.pdf', 'ContentType', 'vector')