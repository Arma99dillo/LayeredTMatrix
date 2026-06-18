function b = MultiLayerSolve(Nlayer,Tmat,a_u,a_d)

%builds the coupled system and solves it using the band gaussian
%elimination algorithm

disp('Building and solving the coupled system')
tic()

%get parameters
M=Tmat{1}.param.M; %they have all the same M
Ndim=2*M+1; 

%get distances between layers and build translation matrices
for j=1:Nlayer-1
    deltaH = Tmat{j}.Hm-Tmat{j+1}.Hp; 

    %get alpha_n and beta_n
    r=1; 
    beta_n_p_2=Tmat{j+1}.beta_n_p; 
    vec=zeros(Ndim,1);
    for n=-M:M
        vec(r)=exp(1i*beta_n_p_2(r)*deltaH);
        r=r+1;
    end

    %build translation matrix
    D{j}=diag(vec);
end

%build coupled system
S=speye(2*Nlayer*Ndim);

%first layer
S(1:Ndim,2*Ndim+1:3*Ndim)=-Tmat{1}.matrix(1:Ndim,1:Ndim)*D{1}; 
S(Ndim+1:2*Ndim,2*Ndim+1:3*Ndim)=-Tmat{1}.matrix(Ndim+1:end,1:Ndim)*D{1}; 

%middle layers
r=2; c=2;
for j=2:Nlayer-1
    S(r*Ndim+1:(r+1)*Ndim,(c-1)*Ndim+1:c*Ndim)=-Tmat{j}.matrix(1:Ndim,Ndim+1:end)*D{j-1};
    S((r+1)*Ndim+1:(r+2)*Ndim,(c-1)*Ndim+1:c*Ndim)=-Tmat{j}.matrix(Ndim+1:end,Ndim+1:end)*D{j-1};
    S(r*Ndim+1:(r+1)*Ndim,(c+2)*Ndim+1:(c+3)*Ndim)=-Tmat{j}.matrix(1:Ndim,1:Ndim)*D{j};
    S((r+1)*Ndim+1:(r+2)*Ndim,(c+2)*Ndim+1:(c+3)*Ndim)=-Tmat{j}.matrix(Ndim+1:end,1:Ndim)*D{j};
    c=c+2; r=r+2;
end

%last layer
S(2*(Nlayer-1)*Ndim+1:2*(Nlayer-1)*Ndim+Ndim,2*(Nlayer-1)*Ndim-Ndim+1:2*(Nlayer-1)*Ndim)=-Tmat{Nlayer}.matrix(1:Ndim,Ndim+1:end)*D{Nlayer-1};
S(2*(Nlayer-1)*Ndim+Ndim+1:end,2*(Nlayer-1)*Ndim-Ndim+1:2*(Nlayer-1)*Ndim)=-Tmat{Nlayer}.matrix(Ndim+1:end,Ndim+1:end)*D{Nlayer-1};

%build rhs
t=zeros(2*Nlayer*Ndim,1);
t(1:Ndim)=Tmat{1}.matrix(1:Ndim,Ndim+1:end)*a_d;
t(Ndim+1:2*Ndim)=Tmat{1}.matrix(Ndim+1:end,Ndim+1:end)*a_d;
t(2*(Nlayer-1)*Ndim+1:2*(Nlayer-1)*Ndim+Ndim)=Tmat{Nlayer}.matrix(1:Ndim,1:Ndim)*a_u;
t(2*(Nlayer-1)*Ndim+Ndim+1:end)=Tmat{Nlayer}.matrix(Ndim+1:end,1:Ndim)*a_u;

%solve system using the Band Gaussian Elimination
b=PentaDiagSolve(S,t,Ndim);

toc()

end