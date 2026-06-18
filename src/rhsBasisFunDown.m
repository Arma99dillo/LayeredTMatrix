function b = rhsBasisFunDown(mesh,param,n)
%computes the rhs for the DtN-TDG method

%get parameters
t=mesh.t; p=mesh.p; K=param.K; delta=param.delta; 
epsilon=param.epsilon; alp=param.alp; LD=mesh.LD; H=param.H; L=param.L; 
nd=param.nd; d=param.d;

%rhs definition
m = size(t,1);
b = zeros(m*nd,1);

k=K*sqrt(epsilon(end)); %wavenumber on the lower boundary
alpha_n=alp+2*pi*n/L;
beta_n=sqrt(k^2-alpha_n^2);

toll=1e-10;

for F=1:size(LD,1) %cycle on the upper boundary
    t1=LD(F,3); %adjiacent element
    p1=p(LD(F,1),:); p2=p(LD(F,2),:); %endpoints
    x1=p1(1); x2=p2(1); %segment endpoints

    for j=1:nd %PW directions cycle
        dj = d(:,j); %PW directions

        %n-th Fourier coefficients of phi_l
        if abs(k*dj(1)-alpha_n) <= toll
            phi_j_n=(1/(L)).*exp(-1i*k*dj(2)*H).*abs(x2-x1);
            int = abs(x2-x1);
        else
            phi_j_n=(1/(L)).*(exp(-1i*k*dj(2)*H)./(1i*(k*dj(1)-alpha_n))).*(exp(1i*(k*dj(1)-alpha_n)*x2)-exp(1i*(k*dj(1)-alpha_n)*x1));
            int = 1./(1i*(alpha_n - k*dj(1))).*(exp(1i*(alpha_n - k*dj(1))*x2)-exp(1i*(alpha_n - k*dj(1))*x1));
        end

        %update rhs
        b((t1-1)*nd+j)=b((t1-1)*nd+j) - 2*1i*beta_n*(1+delta*dj(2))*exp(1i*k*dj(2)*H)*int - (2*1i*delta*L*abs(beta_n)^2*phi_j_n'/k);
    end
end


return