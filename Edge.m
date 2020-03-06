function [E,Eall] = Edge(TR)
    format long;
    points = TR.Points;
    sommets = TR.ConnectivityList;
    E_doublon=zeros(3*size(sommets,1),2);
    for i=1:size(sommets,1)
        E_doublon(3*i-2:3*i,:)=[sommets(i,1),sommets(i,2);sommets(i,2),sommets(i,3);sommets(i,1),sommets(i,3)];
    end
    %E_doublon=[triangles(:,1);triangles(:,2);triangles(:,1)
    %    ,triangles(:,2);triangles(:,3);triangles(:,3)];
    E_doublon_sorted = [min(E_doublon,[],2),max(E_doublon,[],2)];
    Eall=E_doublon;
    %[E,ia1,ic1]=unique(E_doublon);
    E=unique(E_doublon_sorted,'rows','stable');
end