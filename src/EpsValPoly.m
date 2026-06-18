function E = EpsValPoly(mesh,param)
%returns a vector indicating the region where each triangle is located

%get parameters
p=mesh.p; t=mesh.t; vertices=param.vertices;

E=zeros(size(t,1),1);
for j=1:size(t,1)
    bar = (p(t(j,1),:) + p(t(j,2),:) + p(t(j,3),:))/3; %triangle baricenter
    if inpolygon(bar(1),bar(2),vertices(:,1),vertices(:,2))
        E(j)=2;
    else
        E(j)=1;
    end
end

return