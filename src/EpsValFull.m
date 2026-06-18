function E = EpsValFull(mesh,param)
%returns a vector indicating the region where each triangle is located

%get parameters
p=mesh.p; t=mesh.t; v2=param.v2; pos=param.pos;

E=zeros(size(t,1),1);
for j=1:size(t,1)
    bar = (p(t(j,1),:) + p(t(j,2),:) + p(t(j,3),:))/3; %triangle baricenter
    if inpolygon(bar(1),bar(2),v2(:,1),v2(:,2)+pos(2))
        E(j)=2;
    elseif inpolygon(bar(1),bar(2),v2(:,1),v2(:,2)+pos(3))
        E(j)=2;
    elseif inpolygon(bar(1),bar(2),v2(:,1),v2(:,2)+pos(5))
        E(j)=2;
    elseif inpolygon(bar(1),bar(2),v2(:,1),v2(:,2)+pos(8))
        E(j)=2;
    else
        E(j)=1;
    end
end

return