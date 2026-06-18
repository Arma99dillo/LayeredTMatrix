function [mesh,param,u] = FullSolutionVariableShapes(param,LayerArr)

% Solve with DtN-TDG method on full domain

%Global parameters 
param.epsilon=[1, 2, 1]; 

param.pos = LayerArr.pos; % layer center position
param.v1=[2,-0.5;2.5,0.5;3,-0.5]; %triangle
param.v2=[1, 0; 1.5, 0.5; 2, 0; 1.5, -0.5]; %square
param.shape=LayerArr.shape;
param.refl=LayerArr.refl;
param.H=param.pos(1)+1;

mesh.DirObst = true;

%define triangulation
[mesh.p,mesh.t,mesh.I,mesh.B] = MeshVariableShapes(param);

%define matrices for the element edges and divide them in
%internal, upper, lower and left/right edges
[mesh.LI,mesh.LU,mesh.LD,mesh.LR,mesh.LDir] = MeshEdgesDir(mesh,param);

%assign to every mesh element the corresponding epsilon value
mesh.E = EpsValFullVariableShapes(mesh,param);

%assemble and solve linear system
disp('Linear system assembly')
A = MatrixDtNTDG(mesh,param); %system matrix
b = rhsDtNTDG(mesh,param); %system rhs
u=A\b; %solve the system
disp('Solved linear system')


end