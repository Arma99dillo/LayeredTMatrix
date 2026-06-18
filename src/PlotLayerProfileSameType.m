function PlotLayerProfileSameType(LayerArr,Tmat,C)

%plots layer shape and layer limits

%get parameters
%get parameters
L=Tmat{1}.param.L; 
shape=LayerArr.shape; pos=LayerArr.pos; refl=LayerArr.refl;
NLayer=size(shape,1);

fig = gcf; % Get current figure (or use a specific figure handle)
s = numel(findall(fig, 'type', 'axes')); %number of subplots

for l=1:s

    subplot(1,s,l)
    for j=1:NLayer
        plot3([0,L],[Tmat{shape(j)}.Hm+pos(j),Tmat{shape(j)}.Hm+pos(j)],10*ones(1,2),'r--',LineWidth=1.8)
        plot3([0,L],[Tmat{shape(j)}.Hp+pos(j),Tmat{shape(j)}.Hp+pos(j)],10*ones(1,2),'r--',LineWidth=1.8)
        if strcmp(Tmat{shape(j)}.shape.name,'penetrable_poly')
            v1=Tmat{shape(j)}.param.vertices(:,1);
            if refl(j) == 1
                v2=-Tmat{shape(j)}.param.vertices(:,2);
            else
                v2=Tmat{shape(j)}.param.vertices(:,2);
            end
            v2=v2+pos(j);
            v1=[v1;v1(1,:)]; v2=[v2;v2(1,:)];
            plot3(v1,v2,10*ones(size(v1)),'r-',LineWidth=1.8)
        elseif strcmp(Tmat{shape(j)}.shape.name,'dir_poly')
            v1=Tmat{shape(j)}.param.vertices(:,1);
            if refl(j) == 1
                v2=-Tmat{shape(j)}.param.vertices(:,2);
            else
                v2=Tmat{shape(j)}.param.vertices(:,2);
            end
            v2=v2+pos(j);
            v1=[v1;v1(1,:)]; v2=[v2;v2(1,:)];
            plot3(v1,v2,10*ones(size(v1)),'r-',LineWidth=1.8)
        end
        hold on
    end
    axis equal;
    xlim([0,L]); 
    ylim([Tmat{shape(NLayer)}.Hm+pos(NLayer)-C,Tmat{shape(1)}.Hp+pos(1)+C]);

end

end