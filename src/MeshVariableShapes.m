function [p_new,t_new,I,B] = MeshVariableShapes(param)

%Get the parameters
H=param.H; h=param.h; L=param.L; 
pos=param.pos; shape=param.shape; refl=param.refl;

% Create a PDE model
model = createpde();

% Define the geometry 
V = [0,-H; 0,H; L,H; L,-H];  v1=param.v1; v2=param.v2;
n1=size(v1,1); n2=size(v2,1);
g = [3;4;V(:,1);V(:,2)];

geom_list = {g};
for j = 1:size(pos,1)
    if shape(j) == 1
        if refl(j) == 1
            g1 = [2;n1;v1(:,1);-v1(:,2)+pos(j);0;0]; % flipped triangle
        else
            g1 = [2;n1;v1(:,1);v1(:,2)+pos(j);0;0]; % triangle
        end
    else
        if refl(j) == 1
            g1 = [2;n2;v2(:,1);-v2(:,2)+pos(j)]; % flipped square
        else
            g1 = [2;n2;v2(:,1);v2(:,2)+pos(j)]; % square
        end
    end
    geom_list{end+1} = g1; 
end
geom_matrix = [geom_list{:}];

labels = {'g'};
for j = 1:size(pos,1)
    labels{end+1} = sprintf('g%d', j);
end
label_chars = char(labels);

add_terms = {};
sub_terms = {};
for j = 1:size(pos,1)
    if shape(j) == 1
        add_terms{end+1} = sprintf('g%d', j);
    elseif shape(j) == 2
        sub_terms{end+1} = sprintf('g%d', j);
    end
end

% Assemble formula: start with base 'g', then add/subtract groups
set_formula = 'g';
if ~isempty(sub_terms)
    set_formula = [set_formula, '-(' strjoin(sub_terms, '+') ')'];
end
if ~isempty(add_terms)
    set_formula = [set_formula, '+(' strjoin(add_terms, '+') ')'];
end
  
% Add the geometry to the PDE model
geometryFromEdges(model, decsg(geom_matrix, set_formula, label_chars'));


% Generate the mesh with Hmax set to h
mesh = generateMesh(model, 'Hmax', h,'GeometricOrder','linear');


% order triangles and nodes putting the internal nodes 
% before the boundary ones
[p1,e,t1] = meshToPet(mesh); p=p1'; t=t1(1:3,:)';

n_nodes = size(p,1);
n_triangles = size(t,1);


% anti-clockwise orienting of vertices
 for i_t=1:n_triangles
            v=p(t(i_t,:),:);
            if det([v(:,1) v(:,2) ones(3,1)]) < 0
                beep
                t(i_t,[1 2 3])= t(i_t,[2 1 3]);
            end
 end
 
% nodes reordering
old2new = zeros(n_nodes,1);
new2old = zeros(n_nodes,1);
i_int = 1;
i_bound = n_nodes;
for j = 1:n_nodes
    % for every node, I check if it is a boundary or internal node
    if (ismember(j,e(1,(e(6,:)==0) ~= (e(7,:)==0))) || ...
        ismember(j,e(2,(e(6,:)==0) ~= (e(7,:)==0))))
        old2new(j) = i_bound;
        new2old(i_bound) = j;
        i_bound = i_bound - 1;
    else
        old2new(j) = i_int;
        new2old(i_int) = j;
        i_int = i_int + 1;
    end
end

% ordered mesh
p_new = p(new2old,:); 
t_new = old2new(t);
B=(i_bound+1):n_nodes; %boundary indices
I=1:(i_int-1); %internal indices

end