function int = IsInt(T,u,B,t,p,L)
%returns true if in internal element, false otherwise 
toll=1e-6;

if(ismember(t(T,u),B) && ismember(t(T,u+1),B))
    if (abs(p(t(T,u),1)-p(t(T,u+1),1))<=toll) && (abs(p(t(T,u),1))<=toll || abs(p(t(T,u),1)-L)<=toll)
        %left/right boundary
        int=true;
    else
        %boundary element
        int=false;
    end
else
    int=true;
end

return