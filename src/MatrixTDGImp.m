function A = MatrixTDGImp(mesh,param)
%define the system matrix for the TDG method

%get the parameters
t=mesh.t; p=mesh.p; B=mesh.B; K=param.K; LI=mesh.LI; alpha=param.alpha; 
beta=param.beta; delta=param.delta; epsilon=param.epsilon; alp=param.alp; 
E=mesh.E; LU=mesh.LU; LD=mesh.LD; LR=mesh.LR; S=param.L; nd=param.nd; d=param.d;

%sparse matrix allocation
m = size(t,1);
A=spalloc(m*nd,m*nd,10*m*nd);

R =[0,1;-1,0]; %rotation matrix (used for normal vectors)
t = [t t(:,1)]; %modify t for easily deal with the edges

%-----------------------------------
%Internal elements
%-----------------------------------

for T = 1:m 
    %cicle on all the elements to deal with matrix terms associated to
    %basis functions defined on the same element

    k=K*sqrt(epsilon(E(T))); %wavenumber inside the element

    %identify the internal edges
    %include the left and right boundary between the internal edges
    A_aux =zeros(nd,nd); %auxiliary matrix
    for u=1:3
        if IsInt(T,u,B,t,p,S) %check if it's an internal edge
            %cycle on PW directions: this way I include all basis functions
            p1=p(t(T,u),:)'; p2=p(t(T,u+1),:)'; %endpoints
            n = (R*(p2-p1))/norm(p2-p1); %outward normal
            for l = 1:nd
                for j= 1:nd
                    dl = d(:,l); dj = d(:,j); diff=dl-dj; %PW directions and difference
                    A_aux(j,l) = A_aux(j,l) + 1i*k.*(-(1/2).*(dot(dj,n)+dot(dl,n))-beta.*dot(dj,n).*dot(dl,n)-alpha).*phi_int(diff,p1,p2,k);
                end
            end
        end
    end    
    A((T-1)*nd+1:T*nd,(T-1)*nd+1:T*nd)=A((T-1)*nd+1:T*nd,(T-1)*nd+1:T*nd)+A_aux; %update matrix entries
end

%cycle on internal edges to deal with basis functions associated to two
%adjiacent elements
for L=1:size(LI,1)

    t1=LI(L,3); %first element T1
    t2=LI(L,4); %second element T2
    p1=p(LI(L,1),:)'; p2=p(LI(L,2),:)'; %endpoints
    n = (R*(p2-p1))/norm(p2-p1); %outwad normal to T1
    k1=K*sqrt(epsilon(E(t1))); k2=K*sqrt(epsilon(E(t2))); %wavenumbers inside T! and T2
    
    %first cycle: consider T1 for test and T2 for trial   
    A_aux =zeros(nd,nd); %auxiliary matrix
    for l=1:nd %cycle on PW directions
        for j=1:nd
        dl = d(:,l); dj = d(:,j); diff=k1.*dl-k2.*dj; %PW directions and difference
        A_aux(j,l) = A_aux(j,l) + (0.5*(1i*k2*dot(dj,n)+1i*k1*dot(dl,n))+alpha*1i*k2+beta*1i*k1*dot(dl,n)*dot(dj,n)).*phi_int(diff,p1,p2,1);           
        end
    end
    A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)=A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)+A_aux; %update matrix entries
    
    %second cycle: consider T2 for test and T1 for trial 
    %change the normal and the wavenumbers
    n=-n; k1=K*sqrt(epsilon(E(t2))); k2=K*sqrt(epsilon(E(t1)));  
    A_aux =zeros(nd,nd); %auxiliary matrix
    for l=1:nd %cycle on PW directions
        for j=1:nd
        dl = d(:,l); dj = d(:,j); diff=k1.*dl-k2.*dj; %PW directions and difference
        A_aux(j,l) = A_aux(j,l) + (0.5*(1i*k2*dot(dj,n)+1i*k1*dot(dl,n))+alpha*1i*k2+beta*1i*k1*dot(dl,n)*dot(dj,n)).*phi_int(diff,p1,p2,1);           
        end
    end
    A((t1-1)*nd+1:t1*nd,(t2-1)*nd+1:t2*nd)=A((t1-1)*nd+1:t1*nd,(t2-1)*nd+1:t2*nd)+A_aux; %update matrix entries 

end

%-----------------------------------
%Boundary elements
%-----------------------------------

%Upper boundary - Impedance condition
for index=1:size(LU,1)
    t1=LU(index,3); %adjiacent element
    k=K*sqrt(epsilon(E(t1))); %wavenumber in the element
    p1=p(LU(index,1),:)'; p2=p(LU(index,2),:)'; %endpoints
    n = (R*(p2-p1))/norm(p2-p1); %outward normal
    A_aux =zeros(nd,nd); %auxiliary matrix
    for l=1:nd %PW directions cycle
        for j=1:nd
            dl = d(:,l); dj = d(:,j); diff=dl-dj; %PW directions and difference
            A_aux(j,l) = A_aux(j,l) + 1i*k.*((1-delta).*(-1-dot(dj,n))+delta.*dot(dl,n).*(-dot(dj,n)-1)).*phi_int(diff,p1,p2,k);
        end
    end
    A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)=A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)+A_aux; %update matrix
end

%Lower boundary - Impedance condition
for index=1:size(LD,1)
    t1=LD(index,3); %adjiacent element
    k=K*sqrt(epsilon(E(t1))); %wavenumber in the element
    p1=p(LD(index,1),:)'; p2=p(LD(index,2),:)'; %endpoints
    n = (R*(p2-p1))/norm(p2-p1); %outward normal
    A_aux =zeros(nd,nd); %auxiliary matrix
    for l=1:nd %PW directions cycle
        for j=1:nd
            dl = d(:,l); dj = d(:,j); diff=dl-dj; %PW directions and difference
            A_aux(j,l) = A_aux(j,l) + 1i*k.*((1-delta).*(-1-dot(dj,n))+delta.*dot(dl,n).*(-dot(dj,n)-1)).*phi_int(diff,p1,p2,k);
        end
    end
    A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)=A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)+A_aux; %update matrix
end


%Left and right boundary
%deal with all the terms associated to quasi-periodicity: her we consider
%the left and right boundary as the same boundary
c=exp(1i*alp*S); %quasi-periodicity constant
for L=1:size(LR,1)
    t1=LR(L,3); %first element T1 (right)
    t2=LR(L,4); %second element T2 (left)
    p1=p(LR(L,1),:)'; p2=p(LR(L,2),:)'; %endpoints
    n = (R*(p2-p1))/norm(p2-p1); %outward normal to T1
    k=K*sqrt(epsilon(E(t1))); %wavenumber in the elements (it's the same)

    %first cycle: consider T1 for test and T2 for trial   
    A_aux =zeros(nd,nd); %auxiliary matrix
    for l=1:nd %PW directions cycle
        for j=1:nd
        dl = d(:,l); dj = d(:,j); %PW directions
        A_aux(j,l) = A_aux(j,l) - c'.*1i*k.*(-(1/2)*(dot(dj,n)+dot(dl,n))-beta*dot(dj,n)*dot(dl,n)-alpha).*phi_int_bound(dl,dj,p1,p2,k,1,S);            
        end
    end
    A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)=A((t2-1)*nd+1:t2*nd,(t1-1)*nd+1:t1*nd)+A_aux; %update matrix

    n = -n; %change the normal to be outward to T2
    %second cycle: consider T2 for test and T1 for trial  
    A_aux =zeros(nd,nd); %auxiliary matrix
    for l=1:nd %PW directions cycle
        for j=1:nd
        dl = d(:,l); dj = d(:,j); %PW directions
        A_aux(j,l) = A_aux(j,l) - c.*1i*k.*(-(1/2)*(dot(dj,n)+dot(dl,n))-beta*dot(dj,n)*dot(dl,n)-alpha).*phi_int_bound(dl,dj,p1,p2,k,0,S);          
        end
    end
    A((t1-1)*nd+1:t1*nd,(t2-1)*nd+1:t2*nd)=A((t1-1)*nd+1:t1*nd,(t2-1)*nd+1:t2*nd)+A_aux; %update matrix
end

%cycle on Dirichlet boundary if necessary
if mesh.DirObst
    LDir=mesh.LDir;
    for index=1:size(LDir,1)
        t1=LDir(index,3); %adjiacent element
        k=K*sqrt(epsilon(E(t1))); %wavenumber in the element
        p1=p(LDir(index,1),:)'; p2=p(LDir(index,2),:)'; %endpoints
        n = (R*(p2-p1))/norm(p2-p1); %outward normal
        A_aux =zeros(nd,nd); %auxiliary matrix
        for l=1:nd %PW directions cycle
            for j=1:nd
                dl = d(:,l); dj = d(:,j); diff=dl-dj; %PW directions and difference
                A_aux(j,l) = A_aux(j,l) + 1i*k.*(-alpha-dot(dl,n)).*phi_int(diff,p1,p2,k);
            end
        end
        A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)=A((t1-1)*nd+1:t1*nd,(t1-1)*nd+1:t1*nd)+A_aux; %update matrix
    end
end

return