function E = EpsValFullVariableShapes(mesh,param)
%returns a vector indicating the region where each triangle is located

%get parameters
p=mesh.p; t=mesh.t; v1=param.v1; 
pos=param.pos; shape=param.shape; refl=param.refl;

E=zeros(size(t,1),1);
for j=1:size(t,1)
    bar = (p(t(j,1),:) + p(t(j,2),:) + p(t(j,3),:))/3; %triangle baricenter
    cont=0;

    for l=1:size(shape,1)
        if shape(l) == 1
            if refl(l) == 1
                if inpolygon(bar(1),bar(2),v1(:,1),-v1(:,2)+pos(l))
                    cont=1;
                end
            else
                if inpolygon(bar(1),bar(2),v1(:,1),v1(:,2)+pos(l))
                    cont=1;
                end
            end
        end
    end

    if cont == 1
        E(j)=2;
    else
        E(j)=1;
    end
end

return