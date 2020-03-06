function Tri = CreateTriMesh2D(h,lambda,uniform,condition,V)
format long;

wbf = waitbar(0,'Creating 2D uniform mesh');
fd = @(p) sqrt(sum(p.^2,2))-1;
if uniform==0
    fh = @(p) max(5*sqrt(sum(p.^2,2))-1,2);
else
    fh = @(p) ones(size(p,1),1);
end
[p,t] = distmesh(fd,fh,h,[-1,-1;1,1],[],'IT_MAX',100000);

waitbar(.1,wbf,'Performing Delaunay Triangulation');
dt = delaunayTriangulation(lambda*p);

waitbar(.2,wbf,'Retrieving the convex hull');
ch = convexHull(dt);
type(ch)= condition;
x = p(:,1);y=p(:,2);
NbTri = size(dt.ConnectivityList,1);

waitbar(.3,wbf,'Getting all edges of the triangulation');
[E,Eall]=Edge(dt);
NbSeg = size(Eall,1);
Tri={};
Tri{1} = dt.Points;
Tri{2} = dt.ConnectivityList;
Tri{3} = Eall;

waitbar(.35,wbf,'Getting the length of edges');
Long_segment = zeros(1,NbSeg);
for i=1:NbSeg
    Long_segment(1,i) = pdist(dt.Points(Eall(i,:),:),'euclidean');
end
Tri{10} = Long_segment;

waitbar(.4,wbf,'Building the neighborhood of triangles');
TN = 0*dt.ConnectivityList;
for i=1:size(dt.ConnectivityList,1)
    for j=1:size(dt.ConnectivityList,2)
        ID = edgeAttachments(dt,Eall(3*(i-1)+j,1),Eall(3*(i-1)+j,2)); 
        ID=ID{:};
        if size(ID)==1
            TN(i,j)= NaN; 
        elseif ID(1)==i
            TN(i,j) = ID(2);
        else
            TN(i,j) = ID(1);
        end
    end
end

waitbar(.5,wbf,'Adding bordery conditions');
Typesegment = zeros(1,NbSeg);
for k=1:size(ch)-1
    for i=1:NbSeg
        if (ch(k)==Eall(i,1)).*(ch(k+1)==Eall(i,2))
            Typesegment(i)=condition;
        elseif (ch(k+1)==Eall(i,1)).*(ch(k)==Eall(i,2))
            Typesegment(i)=condition;
        end
    end
end
Tri{4} = Typesegment; 

waitbar(.6,wbf,'Finding the circumcenter of the triangles');
Tri{5} = TN;
Centre_tri=circumcenter(dt)';
Tri{6} = Centre_tri;

waitbar(.7,wbf,'Adding distances of triangle with their neighborhood');
Dist_tri = zeros(NbTri,3);
for i=1:NbTri
    for j=1:size(TN,2)
        if isnan(TN(i,j))
            Dist_tri(i,j) = norm(Centre_tri(:,i)'-mean(dt.Points(Eall(3*(i-1)+j,:),:)),2);
        else
            Dist_tri(i,j) = norm((Centre_tri(:,i)-Centre_tri(:,TN(i,j))),2) ;
        end
    end
end
Tri{7} = Dist_tri;

waitbar(.8,wbf,'Creating transmissibility coefficient');
Coef_trans = 0*Dist_tri;
for i=1:NbTri
    Coef_trans(i,1) = Long_segment(1,3*i-2)/Dist_tri(i,1) ;
    Coef_trans(i,2) = Long_segment(1,3*i-1)/Dist_tri(i,2) ;
    Coef_trans(i,3) = Long_segment(1,3*i)/Dist_tri(i,3) ;
end
Tri{8} = Coef_trans;

waitbar(.85,wbf,'Adding surface to each triangle');
Volume = zeros(1,NbTri);
for i=1:NbTri
    s = 0.5*sum(Long_segment(3*i-2:3*i));
    Volume(1,i) = sqrt(s*(s-Long_segment(3*i-2))*(s-Long_segment(3*i-1))*(s-Long_segment(3*i)));
end
Tri{9} = Volume;

waitbar(.9,wbf,'Approximating convection-term');
Vites_seg = zeros(NbSeg,2);
Normales_seg = zeros(NbSeg,2);
Prodvitnor_seg = zeros(NbSeg,1);
for i=1:NbSeg
    Vites_seg(i,:) = V(0.5*(dt.Points(Eall(i,1),:)+dt.Points(Eall(i,2),:)));
    Normales_seg(i,:) = (0.5*(dt.Points(Eall(i,1),:)+dt.Points(Eall(i,2),:))-Centre_tri(:,floor((i-1)/3)+1)')/(norm((0.5*(dt.Points(Eall(i,1),:)+dt.Points(Eall(i,2),:))-Centre_tri(:,floor((i-1)/3)+1)'),2));
    Prodvitnor_seg(i,1) = -Vites_seg(i,:)*Normales_seg(i,:)';
end
Tri{11} = Vites_seg;
Tri{12} = Normales_seg;
Tri{13} = Prodvitnor_seg;

waitbar(1,wbf,'Done');
close(wbf);
end