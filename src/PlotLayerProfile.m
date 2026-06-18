function PlotLayerProfile(Tmat,C)

%plots layer shape and layer limits

fig = gcf; % Get current figure (or use a specific figure handle)
s = numel(findall(fig, 'type', 'axes')); %number of subplots

L=Tmat{1}.param.L;

for j=1:size(Tmat,2)

    TmatTemp=Tmat{j};
    if strcmp(TmatTemp.param.name,'u_shape')
        d=(TmatTemp.Hp+TmatTemp.Hm)/2-1;

        for l=1:s
            subplot(1,s,l)
            plot3([0,L/4],[d,d],10*ones(1,2),'r-',LineWidth=1.8)
            plot3([3*L/4,L],[d,d],10*ones(1,2),'r-',LineWidth=1.8)
            plot3([L/4,L/4],[d,d+2],10*ones(1,2),'r-',LineWidth=1.8)
            plot3([3*L/4,3*L/4],[d,d+2],10*ones(1,2),'r-',LineWidth=1.8)
            plot3([L/4,3*L/4],[d+2,d+2],10*ones(1,2),'r-',LineWidth=1.8)
            plot3([0,L],[TmatTemp.Hm,TmatTemp.Hm],10*ones(1,2),'r--',LineWidth=1.8)
            plot3([0,L],[TmatTemp.Hp,TmatTemp.Hp],10*ones(1,2),'r--',LineWidth=1.8)
            hold on;
        end

    elseif strcmp(TmatTemp.param.name,'penetrable_poly') || strcmp(TmatTemp.param.name,'dir_poly')
        v1=TmatTemp.param.vertices(:,1); v2=TmatTemp.param.vertices(:,2)+(TmatTemp.Hp-TmatTemp.param.H);
        v1=[v1;v1(1,:)]; v2=[v2;v2(1,:)];

        for l=1:s
            subplot(1,s,l)
            plot3(v1,v2,10*ones(size(v1)),'r-',LineWidth=1.8)
            plot3([0,L],[TmatTemp.Hp,TmatTemp.Hp],10*ones(1,2),'r--',LineWidth=1.8)
            plot3([0,L],[TmatTemp.Hm,TmatTemp.Hm],10*ones(1,2),'r--',LineWidth=1.8)
            hold on;
        end

    elseif strcmp(TmatTemp.param.name,'double_rectangle')
        d=(TmatTemp.Hp+TmatTemp.Hm)/2;

        for l=1:s
            subplot(1,s,l)
            plot3([0,L],[d,d],10*ones(1,2),'r-',LineWidth=1.8)
            plot3([0,L],[TmatTemp.Hp,TmatTemp.Hp],10*ones(1,2),'r--',LineWidth=1.8)
            plot3([0,L],[TmatTemp.Hm,TmatTemp.Hm],10*ones(1,2),'r--',LineWidth=1.8)
            hold on;
        end

    end

end

for l=1:s
    subplot(1,s,l)
    xlim([0,L]); ylim([Tmat{end}.Hm-C,Tmat{1}.Hp+C]); view(2)
    hold on
end

end