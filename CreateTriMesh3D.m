function Tetra = CreateTriMesh3D(h,lambda,uniform,condition,V)
wbf = waitbar(0,'Creating 3D uniform mesh');
fd = @(p) sqrt(sum(p.^2,2)) - 1;
if uniform==0
    fh = @(p) max(2*sqrt(sum(p.^2,2))-1,2);%ones(size(p,1),1);
elseif uniform==1
    fh = @(p) ones(size(p,1),1);
end
[p,t] = distmesh( fd, fh, h, [-1,-1,-1;1,1,1],[],'IT_MAX',100000 );

waitbar(.2,wbf,'Performing Delaunay Triangulation');
DT = delaunayTriangulation(lambda*p);
Tetra = {};
Tetra{1} = DT.Points;
Tetra{2} = DT.ConnectivityList;

waitbar(.3,wbf,'Creating the list of the triangles');
sommets = DT.ConnectivityList;
Triang=zeros(4*size(sommets,1),3);
for i=1:size(sommets,1)
    Triang(4*i-3:4*i,:)=[sommets(i,1),sommets(i,2),sommets(i,3);sommets(i,2),sommets(i,3),sommets(i,4);sommets(i,1),sommets(i,3),sommets(i,4);sommets(i,1),sommets(i,2),sommets(i,4)];
end
Tetra{3} = Triang;

waitbar(.4,wbf,'Getting the neighborhood of tetrahedron');
TNbad = neighbors(DT);
TN=0*TNbad;
for i=1:size(TNbad,1)
    for k=1:4
        if isnan(TNbad(i,1))
        else
            if (size(intersect(Triang(4*(i-1)+k,:),DT.ConnectivityList(TNbad(i,1),:)),2)==3)
                TN(i,k) = TNbad(i,1);
            end
        end
        if isnan(TNbad(i,2))
        else
            if (size(intersect(Triang(4*(i-1)+k,:),DT.ConnectivityList(TNbad(i,2),:)),2)==3)
                TN(i,k) = TNbad(i,2);
            end
        end
        if isnan(TNbad(i,3))
        else
            if (size(intersect(Triang(4*(i-1)+k,:),DT.ConnectivityList(TNbad(i,3),:)),2)==3)
                TN(i,k) = TNbad(i,3);
            end
        end
        if isnan(TNbad(i,4))
        else
            if (size(intersect(Triang(4*(i-1)+k,:),DT.ConnectivityList(TNbad(i,4),:)),2)==3)
                TN(i,k) = TNbad(i,4);
            end
        end
        if (TN(i,k)==0)
            TN(i,k) = NaN;
        end
    end
end
Tetra{4} = TN;

waitbar(.5,wbf,'Giving bordery conditions');
[K,v] = convexHull(DT);
Typetriang = zeros(size(Triang,1),1);
for i=1:size(Triang,1)
    for k=1:size(K,1)
        if (min(Triang(i,:))==K(k,1))
            if size(intersect(Triang(i,:),K(k,:)),2)==3
                Typetriang(i,1) = condition;
            end
        end
    end
end
Tetra{5} = Typetriang;

waitbar(.6,wbf,'Getting the surface of triangles');
Surftriang = zeros(size(Triang,1),1);
for i=1:size(Triang,1)
    a = pdist(DT.Points(Triang(i,1:2),:),'euclidean');
    b = pdist(DT.Points(Triang(i,2:3),:),'euclidean');
    c = pdist(DT.Points(Triang(i,[1,3]),:),'euclidean');
    s = 0.5*(a+b+c);
    Surftriang(i,1) = sqrt(s*(s-a)*(s-b)*(s-c));
end
Tetra{6} = Surftriang;

waitbar(.65,wbf,'Getting the circumcenter of tetrahedron');
Centre_tetra=circumcenter(DT);
Tetra{7} = Centre_tetra;

waitbar(.7,wbf,'Adding the distance of tetrahedrons with their neighborhood');
Dist_tetra = 0*TN;
for i=1:size(TN,1)
    for j=1:size(TN,2)
        if isnan(TN(i,j))
            Dist_tetra(i,j) = pdist([Centre_tetra(i,:);center(DT.Points(Triang(4*(i-1)+j,1),:),DT.Points(Triang(4*(i-1)+j,2),:),DT.Points(Triang(4*(i-1)+j,3),:),'circumcenter')'],'euclidean');%mean(DT.Points(Triang(4*(i-1)+j,:),:))
        else 
            Dist_tetra(i,j) = pdist([Centre_tetra(i,:);Centre_tetra(TN(i,j),:)],'euclidean');
        end
    end
end
Tetra{8} = Dist_tetra;

waitbar(.75,wbf,'Creating the transmissibility coefficient');
Coef_trans = 0*Dist_tetra;
for i=1:size(Coef_trans,1)
   for j=1:size(Coef_trans,2)
       Coef_trans(i,j) = Surftriang(4*(i-1)+j,1)/Dist_tetra(i,j);
   end
end
Tetra{9} = Coef_trans;

waitbar(.8,wbf,'Getting the volume of tetrahedrons');
Volume = zeros(size(Centre_tetra,1),1);
for i=1:size(Centre_tetra,1)
    d12 = pdist(DT.Points(DT.ConnectivityList(i,1:2),:));
    d23 = pdist(DT.Points(DT.ConnectivityList(i,2:3),:));
    d34 = pdist(DT.Points(DT.ConnectivityList(i,3:4),:));
    d13 = pdist(DT.Points(DT.ConnectivityList(i,[1,3]),:));
    d14 = pdist(DT.Points(DT.ConnectivityList(i,[1,4]),:));
    d24 = pdist(DT.Points(DT.ConnectivityList(i,[2,4]),:));
    Mat = [0,d12^2,d13^2,d14^2,1;
        d12^2,0,d23^2,d24^2,1;
        d13^2,d23^2,0,d34^2,1;
        d14^2,d24^2,d34^2,0,1;
        1,1,1,1,0];
    Volume(i,1) = sqrt(1/288*det(Mat));
end
Tetra{10} = Volume;

waitbar(.9,wbf,'Approximating convective properties');
Vites_surf = zeros(size(Triang,1),3);
Normales_surf = zeros(size(Triang,1),3);
Prodvitnor_surf = zeros(size(Triang,1),1);
for i=1:size(Triang,1)
    Vites_surf(i,:) = V(center(DT.Points(Triang(i,1),:),DT.Points(Triang(i,2),:),DT.Points(Triang(i,3),:),'circumcenter')');
    Normales_surf(i,:) = (center(DT.Points(Triang(i,1),:),DT.Points(Triang(i,2),:),DT.Points(Triang(i,3),:),'circumcenter')'-Centre_tetra(floor((i-1)/4)+1,:))/(norm((center(DT.Points(Triang(i,1),:),DT.Points(Triang(i,2),:),DT.Points(Triang(i,3),:),'circumcenter')'-Centre_tetra(floor((i-1)/4)+1,:)),2));
    Prodvitnor_surf(i,1) = -Vites_surf(i,:)*Normales_surf(i,:)';
end
Tetra{11} = Vites_surf;
Tetra{12} = Normales_surf;
Tetra{13} = Prodvitnor_surf;

waitbar(1,wbf,'Done');
close(wbf);
end