function [Tmat,s] = MultiTestVariableShape(param,LayerArr)

% Define layer shapes
NShape=2; LayerShape=cell(NShape,1);

%First layer
LayerShape{1}.Hp=1; LayerShape{1}.Hm=-1; LayerShape{1}.epsilon=[1, 2, 1]; 
LayerShape{1}.name='penetrable_poly'; LayerShape{1}.vertices=[2,-0.5;2.5,0.5;3,-0.5];

%Second layer
LayerShape{2}.Hp=1; LayerShape{2}.Hm=-1; LayerShape{2}.epsilon=1; 
LayerShape{2}.name='dir_poly'; LayerShape{2}.vertices=[1, 0; 1.5, 0.5; 2, 0; 1.5, -0.5];

% Mesh and matrices
parfor j=1:NShape
    [Tmat{j}.mesh,Tmat{j}.param] = GenerateMesh(param,LayerShape{j});
    Tmat{j}.Hp=LayerShape{j}.Hp; Tmat{j}.Hm=LayerShape{j}.Hm;
    Tmat{j}.shape=LayerShape{j};
end

parfor j=1:NShape
    disp(['T-Matrix ', num2str(j), ' computation'])
    Tmat{j} = BuildTMatrix(Tmat{j});
end

% Build coupled linear system
M=Tmat{1}.param.M; Ndim=2*M+1;

nsol=0;
a_d=zeros(Ndim,1); a_d(nsol+M+1)=exp(-1i*(Tmat{1}.Hp+LayerArr.pos(1))*Tmat{1}.beta_n_p(nsol+M+1));
a_u=zeros(Ndim,1);

s = MultiLayerSolveSameType(LayerArr,Tmat,a_u,a_d);

end