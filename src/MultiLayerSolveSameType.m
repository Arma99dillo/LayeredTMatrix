function b = MultiLayerSolveSameType(LayerArr,Tmat,a_u,a_d)

%builds the coupled system and solves it using the band gaussian
%elimination algorithm when there are repeated shapes

disp('Building and solving the coupled system')
tic()

%get parameters
shape=LayerArr.shape; pos=LayerArr.pos; refl=LayerArr.refl;
NLayer=size(shape,1);
M=Tmat{1}.param.M; %they have all the same M
Ndim=2*M+1; 

%get distances between layers and build translation matrices
for j=1:NLayer-1
    deltaH = (Tmat{shape(j)}.Hm+pos(j))-(Tmat{shape(j+1)}.Hp+pos(j+1)); 

    %get alpha_n and beta_n
    r=1; 
    beta_n_p_2=Tmat{shape(j+1)}.beta_n_p; 
    vec=zeros(Ndim,1);
    for n=-M:M
        vec(r)=exp(1i*beta_n_p_2(r)*deltaH);
        r=r+1;
    end

    %build translation matrix
    D{j}=diag(vec);
end

%build coupled system
S=speye(2*NLayer*Ndim); t=zeros(2*NLayer*Ndim,1);

%reflect T-matrix if necessary
if refl(1) == 1
    T = ReflectTmat(Tmat{shape(1)}).matrix;
else
    T = Tmat{shape(1)}.matrix;
end

%first layer
S(1:Ndim,2*Ndim+1:3*Ndim)=-T(1:Ndim,1:Ndim)*D{1}; 
S(Ndim+1:2*Ndim,2*Ndim+1:3*Ndim)=-T(Ndim+1:end,1:Ndim)*D{1}; 

t(1:Ndim)=T(1:Ndim,Ndim+1:end)*a_d;
t(Ndim+1:2*Ndim)=T(Ndim+1:end,Ndim+1:end)*a_d;

%middle layers
r=2; c=2;
for j=2:NLayer-1

    if refl(j) == 1
        T = ReflectTmat(Tmat{shape(j)}).matrix;
    else
        T = Tmat{shape(j)}.matrix;
    end
    S(r*Ndim+1:(r+1)*Ndim,(c-1)*Ndim+1:c*Ndim)=-T(1:Ndim,Ndim+1:end)*D{j-1};
    S((r+1)*Ndim+1:(r+2)*Ndim,(c-1)*Ndim+1:c*Ndim)=-T(Ndim+1:end,Ndim+1:end)*D{j-1};
    S(r*Ndim+1:(r+1)*Ndim,(c+2)*Ndim+1:(c+3)*Ndim)=-T(1:Ndim,1:Ndim)*D{j};
    S((r+1)*Ndim+1:(r+2)*Ndim,(c+2)*Ndim+1:(c+3)*Ndim)=-T(Ndim+1:end,1:Ndim)*D{j};
    c=c+2; r=r+2;
end
if refl(NLayer) == 1
    T = ReflectTmat(Tmat{shape(NLayer)}).matrix;
else
    T = Tmat{shape(NLayer)}.matrix;
end

%last layer
S(2*(NLayer-1)*Ndim+1:2*(NLayer-1)*Ndim+Ndim,2*(NLayer-1)*Ndim-Ndim+1:2*(NLayer-1)*Ndim)=-T(1:Ndim,Ndim+1:end)*D{NLayer-1};
S(2*(NLayer-1)*Ndim+Ndim+1:end,2*(NLayer-1)*Ndim-Ndim+1:2*(NLayer-1)*Ndim)=-T(Ndim+1:end,Ndim+1:end)*D{NLayer-1};

t(2*(NLayer-1)*Ndim+1:2*(NLayer-1)*Ndim+Ndim)=T(1:Ndim,1:Ndim)*a_u;
t(2*(NLayer-1)*Ndim+Ndim+1:end)=T(Ndim+1:end,1:Ndim)*a_u;

%Solve system with Band Gaussian Elimination
b=PentaDiagSolve(S,t,Ndim);

toc()

end