function rhs = rhsTDGImp(mesh,Tmat,b_u_p,b_d_p,b_u_m,b_d_m,deltaH_p,deltaH_m)
%computes the rhs for the TDG method

%get parameters
param=Tmat.param; Hp=Tmat.Hp; Hm=Tmat.Hm;
t=mesh.t; p=mesh.p; K=param.K; delta=param.delta; 
epsilon=param.epsilon; LU=mesh.LU; LD=mesh.LD;
nd=param.nd; d=param.d; M=param.M;
alpha_n=Tmat.alpha_n; beta_n_p=Tmat.beta_n_p; beta_n_m=Tmat.beta_n_m;

%rhs definition
rhs = zeros(size(t,1)*nd,1);

toll=1e-10;

k=K*sqrt(epsilon(1)); normal=[0;1]; %wavenumber and normal on the upper boundary
for F=1:size(LU,1) %cycle on the upper boundary
    t1=LU(F,3); %adjiacent element
    p1=p(LU(F,1),:)'; p2=p(LU(F,2),:)'; %endpoints
    x1=p2(1); x2=p1(1); %segment endpoints

    for j=1:nd %PW directions cycle
        dj = d(:,j); %PW direction
        sum=0;

        %cycle on N to compute the integral
        v=1;
        for n=-M:M
            if abs(k*dj(1)-alpha_n(v)) <= toll
                int = abs(x2-x1);
            else
                int = 1./(1i*(alpha_n(v) - k*dj(1))).*(exp(1i*(alpha_n(v) - k*dj(1))*x2)-exp(1i*(alpha_n(v) - k*dj(1))*x1));
            end
            sum=sum+(b_u_p(v)*(1i*beta_n_p(v)-1i*k)+b_d_p(v)*(-1i*beta_n_p(v)-1i*k)*exp(1i*beta_n_p(v)*deltaH_p))*int;
            v=v+1;
        end

        %update rhs
        rhs((t1-1)*nd+j)=rhs((t1-1)*nd+j) + (delta*(-dot(dj,normal)-1)+1)*exp(-1i*k*Hp*dj(2))*sum;
    end
end


k=K*sqrt(epsilon(end)); normal=[0;-1]; %wavenumber and normal on the lower boundary
for F=1:size(LD,1) %cycle on the upper boundary
    t1=LD(F,3); %adjiacent element
    p1=p(LD(F,1),:)'; p2=p(LD(F,2),:)'; %endpoints
    x1=p1(1); x2=p2(1); %segment endpoints

    for j=1:nd %PW directions cycle
        dj = d(:,j); %PW direction
        sum=0;

        %cycle on N to compute the integral
        v=1;
        for n=-M:M
            if abs(k*dj(1)-alpha_n(v)) <= toll
                int = abs(x2-x1);
            else
                int = 1./(1i*(alpha_n(v) - k*dj(1))).*(exp(1i*(alpha_n(v) - k*dj(1))*x2)-exp(1i*(alpha_n(v) - k*dj(1))*x1));
            end
            sum=sum+(b_d_m(v)*(1i*beta_n_m(v)-1i*k)+b_u_m(v)*(-1i*beta_n_m(v)-1i*k)*exp(1i*beta_n_m(v)*deltaH_m))*int;
            v=v+1;
        end

        %update rhs
        rhs((t1-1)*nd+j)=rhs((t1-1)*nd+j) + (delta*(-dot(dj,normal)-1)+1)*exp(-1i*k*Hm*dj(2))*sum;
    end
end


return