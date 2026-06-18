function [mesh,param] = GenerateMesh(param, domain)
%Define the triangulation on the domain specified 

if strcmp(domain.name,'double_rectangle')
    param.H=(domain.Hp-domain.Hm)/2;
    param.epsilon=domain.epsilon;
    param.name=domain.name;
    mesh.DirObst = false; %no Dirichlet obstacle

    %define triangulation - MUST BE PERIODIC
    [ mesh.p,mesh.t,mesh.I,mesh.B ] = MeshDouble(param);

    %define matrices for the element edges and divide them in
    %internal, upper, lower and left/right edges
    [mesh.LI,mesh.LU,mesh.LD,mesh.LR] = MeshEdges(mesh,param);

    %assign to every mesh element the corresponding epsilon value
    mesh.E = EpsValDouble(mesh);

elseif strcmp(domain.name,'u_shape')
    param.H=(domain.Hp-domain.Hm)/2; 
    param.epsilon=domain.epsilon;
    param.name=domain.name;
    mesh.DirObst = false; %no Dirichlet obstacle

    %define triangulation - MUST BE PERIODIC
    [ mesh.p,mesh.t,mesh.I,mesh.B ] = MeshUShape(param);

    %define matrices for the element edges and divide them in
    %internal, upper, lower and left/right edges
    [mesh.LI,mesh.LU,mesh.LD,mesh.LR] = MeshEdges(mesh,param);

    %assign to every mesh element the corresponding epsilon value
    mesh.E = EpsValUShape(mesh,param);

elseif strcmp(domain.name,'dir_poly')
    param.H=(domain.Hp-domain.Hm)/2;
    param.epsilon=domain.epsilon; param.vertices=domain.vertices;
    param.name=domain.name;
    mesh.DirObst = true; %there is a Dirichlet obstacle

    %define triangulation - MUST BE PERIODIC
    [ mesh.p,mesh.t,mesh.I,mesh.B ] = MeshDirPoly(param);

    %define matrices for the element edges and divide them in
    %internal, upper, lower and left/right edges
    [mesh.LI,mesh.LU,mesh.LD,mesh.LR,mesh.LDir] = MeshEdgesDir(mesh,param);

    %assign to every mesh element the corresponding epsilon value
    mesh.E = ones(size(mesh.t,1),1);    

elseif strcmp(domain.name,'penetrable_poly')
    param.H=(domain.Hp-domain.Hm)/2; 
    param.epsilon=domain.epsilon; param.vertices=domain.vertices;
    param.name=domain.name;
    mesh.DirObst = false; %no Dirichlet obstacle

    %define triangulation - MUST BE PERIODIC
    [ mesh.p,mesh.t,mesh.I,mesh.B ] = MeshPenetrablePoly(param);

    %define matrices for the element edges and divide them in
    %internal, upper, lower and left/right edges
    [mesh.LI,mesh.LU,mesh.LD,mesh.LR] = MeshEdges(mesh,param);

    %assign to every mesh element the corresponding epsilon value
    mesh.E = EpsValPoly(mesh,param);

else
    error('Wrong domain!')
end

return