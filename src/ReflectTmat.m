function Tmat_refl = ReflectTmat(Tmat)

%reflects the mesh and T-matrix

Tmat_refl = Tmat;

% Reflect across y = 0
Tmat_refl.mesh.p(:,2) = -Tmat_refl.mesh.p(:,2);

% Reverse triangle orientation so normals stay consistent
Tmat_refl.mesh.t(:,[2 3]) = Tmat_refl.mesh.t(:,[3 2]);

%Flip epsilon and E
Tmat_refl.param.epsilon = fliplr(Tmat_refl.param.epsilon);
Tmat_refl.mesh.E = size(Tmat_refl.param.epsilon,2) + 1 - Tmat_refl.mesh.E;

%Change left, right, upper and lower boundary
[Tmat_refl.mesh.LI,Tmat_refl.mesh.LU,Tmat_refl.mesh.LD,Tmat_refl.mesh.LR] = MeshEdges(Tmat_refl.mesh,Tmat_refl.param);

%Change T-matrix
Ndim = 2*Tmat.param.M+1;
R = [zeros(Ndim,Ndim), eye(Ndim); eye(Ndim), zeros(Ndim,Ndim),];
Tmat_refl.matrix = R*Tmat_refl.matrix*R;
temp = Tmat_refl.beta_n_p;
Tmat_refl.beta_n_p = Tmat_refl.beta_n_m;  Tmat_refl.beta_n_m = temp;

end