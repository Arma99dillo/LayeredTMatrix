function Tmat = BuildTMatrix(Tmat)
%builds the T-matrix using basis functions Psi_n as incident fields

%get parameters
param=Tmat.param; mesh=Tmat.mesh;
M=param.M; LU=mesh.LU; LD=mesh.LD; nd=param.nd;

disp('T-Matrix computation')
tic()

%build TDG stiffness matrix 
[A,FCp,FCm,Tmat.alpha_n,Tmat.beta_n_p,Tmat.beta_n_m] = MatrixDtNTDG(mesh,param); 

%build TDG rhs for every Psi_n
b=zeros(size(A,1), 4*M+2); j=1;
for n=-M:M
    b(:,j)=rhsBasisFunUp(mesh,param,n);
    b(:,j+1)=rhsBasisFunDown(mesh,param,n);
    j=j+2;
end

u_TOT=A\b; %solve linear systems

%build T-matrix
T = zeros(4*M+2,4*M+2);
l=1; v=1;
for n=-M:M
    b_up = zeros(2*M+1,1); b_down = zeros(2*M+1,1);
    u=u_TOT(:,v);
    for F=1:size(LU,1) 
        t=LU(F,3); %corresponding triangle
        u_coeff=u((t-1)*nd+1:t*nd); %solution coefficients
        for j=1:nd
            FC_coeff=FCp((2*M+1)*nd*(F-1)+(j-1)*(2*M+1)+1 : (2*M+1)*nd*(F-1)+j*(2*M+1)); %Fourier coefficients
            b_up=b_up+u_coeff(j)*FC_coeff;
        end
    end
    b_up(n+M+1)=b_up(n+M+1)-1;
    for F=1:size(LD,1) 
        t=LD(F,3); %corresponding triangle
        u_coeff=u((t-1)*nd+1:t*nd); %solution coefficients
        for j=1:nd
            FC_coeff=FCm((2*M+1)*nd*(F-1)+(j-1)*(2*M+1)+1 : (2*M+1)*nd*(F-1)+j*(2*M+1)); %Fourier coefficients
            b_down=b_down+u_coeff(j)*FC_coeff;
        end
    end
    T(:,2*M+1+l)=[b_up; b_down];
    v=v+1;

    b_up = zeros(2*M+1,1); b_down = zeros(2*M+1,1);
    u=u_TOT(:,v);
    for F=1:size(LU,1) 
        t=LU(F,3); %corresponding triangle
        u_coeff=u((t-1)*nd+1:t*nd); %solution coefficients
        for j=1:nd
            FC_coeff=FCp((2*M+1)*nd*(F-1)+(j-1)*(2*M+1)+1 : (2*M+1)*nd*(F-1)+j*(2*M+1)); %Fourier coefficients
            b_up=b_up+u_coeff(j)*FC_coeff;
        end
    end
    for F=1:size(LD,1) 
        t=LD(F,3); %corresponding triangle
        u_coeff=u((t-1)*nd+1:t*nd); %solution coefficients
        for j=1:nd
            FC_coeff=FCm((2*M+1)*nd*(F-1)+(j-1)*(2*M+1)+1 : (2*M+1)*nd*(F-1)+j*(2*M+1)); %Fourier coefficients
            b_down=b_down+u_coeff(j)*FC_coeff;
        end
    end
    b_down(n+M+1)=b_down(n+M+1)-1;
    T(:,l)=[b_up; b_down];
    l=l+1; v=v+1;
end

Tmat.matrix=T;

toc()

return