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

%% Cycle on number of layers
ni = 3; nf = 45; %min and max number of layers
t = zeros((nf-ni)/2+1,2);

v=1; 
for l=ni:2:nf

    disp(['Solving with ', num2str(l), ' obstacles'])

    LayerArr.pos = (3*v:-3:-3*v)'; % layer center position
    LayerArr.shape = randi([1 2], l, 1); % layer shape
    LayerArr.refl = randi([0 1], l, 1); %layer reflection
    
    %solve with "monolithic" solver
    disp('DtN-TDG method on full domain')
    tic
    [TDGmesh,TDGparam,u] = FullSolutionVariableShapes(param,LayerArr);
    t(v,1) = toc;


    %solve with T-matrix
    disp('T-matrix method')
    tic
    [Tmat,s] = MultiTestVariableShape(param,LayerArr);
    t(v,2) = toc;

    v=v+1;
    
end

%% Convergence plot
figure()
semilogy(ni:2:nf,t(:,1),'o-',ni:2:nf,t(:,2),'*-','LineWidth',1.2); grid
xlim([ni nf]); xticks(5:5:nf); 
LL = legend('Full TDG','T-matrix','FontSize', 14);
set(LL, 'Interpreter', 'latex');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14,'TickLabelInterpreter', 'latex')
xlabel('Number of layers','FontSize',18, 'Interpreter','latex')
ylabel('Computational time','FontSize',18, 'Interpreter','latex')

%% Save plot
exportgraphics(gcf, 'figs/CompTime.pdf', 'ContentType', 'vector')