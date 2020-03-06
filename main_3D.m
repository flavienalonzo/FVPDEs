h = 0.25;lambda=1;uniform=1;condition=2;
V = @(x) x*[0,-1,0;1,0,0;0,0,1];

%% 
Tetra = CreateTriMesh3D(h,lambda,uniform,condition,V);
save(['Tetra_3D_',num2str(h),'_',num2str(lambda),'_',num2str(uniform),'_',num2str(condition),'_',num2str(length(func2str(V))),'.mat'],'Tetra');

%% 
load(['Tetra_3D_',num2str(h),'_',num2str(lambda),'_',num2str(uniform),'_',num2str(condition),'_',num2str(length(func2str(V))),'.mat']);
equation=5;

%% 
Points = Tetra{1};
ConnectivityList = Tetra{2};
Triang = Tetra{3};
TN = Tetra{4};
Typetriang = Tetra{5};
Surftriang = Tetra{6};
Centre_tetra = Tetra{7};
Dist_tetra = Tetra{8};
Coef_trans = Tetra{9};
Volume = Tetra{10};

%% 
u0 = zeros(size(Centre_tetra,1),1);
for i=1:size(Centre_tetra,1)
    if pdist([Centre_tetra(i,:);lambda*[0.25,0.25,0]])<2*h
        u0(i,1) = 1;
    end
end
Delta_t = 0.001;

%% Création d'un maillage
h = 0.15;
fd = @(p) sqrt(sum(p.^2,2)) - 1;
fh = @(p) max(2*sqrt(sum(p.^2,2))-1,2);%ones(size(p,1),1);
[p,t] = distmesh( fd, fh, h, [-1,-1,-1;1,1,1],[],'IT_MAX',100000 );
f = [t(:,[1:3]); t(:,[1,2,4]); t(:,[2,3,4]); t(:,[3,1,4])];
patch( 'vertices', p, 'faces', f, 'facecolor', [.9, .9, .9] );
lambda=5;
DT = delaunayTriangulation(lambda*p);
%figure;
%tetramesh(DT,'FaceAlpha',0.3);
Tetra = {};
Tetra{1} = DT.Points;
Tetra{2} = DT.ConnectivityList;

%% Création de la liste des segments du maillage
E = edges(DT);

%% Création de la liste des triangles du maillage
sommets = DT.ConnectivityList;
Triang=zeros(4*size(sommets,1),3);
for i=1:size(sommets,1)
    Triang(4*i-3:4*i,:)=[sommets(i,1),sommets(i,2),sommets(i,3);sommets(i,2),sommets(i,3),sommets(i,4);sommets(i,1),sommets(i,3),sommets(i,4);sommets(i,1),sommets(i,2),sommets(i,4)];
end
Tetra{3} = Triang;
%% Listes des tétraèdres voisins
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

%% Ajout du type de triangle 
[K,v] = convexHull(DT);
condition = input('Quel condition au bord ? 1) Dirichlet 2) Neumann ');
%trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3))
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
if (condition==2)
    Typetriang(1:size(Triang,1),1) = 0.5*Typetriang(1:size(Triang,1),1);
end
Tetra{5} = Typetriang;

%% Ajout de la surface des triangles
Surftriang = zeros(size(Triang,1),1);
for i=1:size(Triang,1)
    a = pdist(DT.Points(Triang(i,1:2),:),'euclidean');
    b = pdist(DT.Points(Triang(i,2:3),:),'euclidean');
    c = pdist(DT.Points(Triang(i,[1,3]),:),'euclidean');
    s = 0.5*(a+b+c);
    Surftriang(i,1) = sqrt(s*(s-a)*(s-b)*(s-c));
end
Tetra{6} = Surftriang;

%% Ajout des centres de chaque tétraèdre
Centre_tetra=circumcenter(DT);
% figure;
% plot3(Centre_tetra(:,1),Centre_tetra(:,2),Centre_tetra(:,3),'+');
Tetra{7} = Centre_tetra;

%% Ajout des distances entre tétraèdre
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

%% Création coefficient de transmissibilité
Coef_trans = 0*Dist_tetra;
for i=1:size(Coef_trans,1)
   for j=1:size(Coef_trans,2)
       Coef_trans(i,j) = Surftriang(4*(i-1)+j,1)/Dist_tetra(i,j);
   end
end
Tetra{9} = Coef_trans;

%% Ajout volume des tétraèdres
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

%% Introduction au problème
equation = input('Quelle équation ? 1) Diffusion 2) Convection 3) Diffusion-Convection ');
Delta_t=0;
if equation==2||equation==3
    vitesse = input('Quelle vitesse de convection ? 1) V=Rz(90°)*x 2) V=Rx(90°).Rz(90°)*x 3) V=Rx(90°).Ry(90°).Rz(90°)*x/norm(x,2) ');
    if vitesse==1
        V = @(x) x*[0,-1,0;1,0,0;0,0,1];
    elseif vitesse==2
        V = @(x) x*[0,-1,0;1,0,0;0,0,1]*[1,0,0;0,0,-1;0,1,0];
    elseif vitesse==3
        V = @(x) x*[0,-1,0;1,0,0;0,0,1]*[0,0,1;0,1,0;-1,0,0]*[1,0,0;0,0,-1;0,1,0]/norm(x,2);
    end
    cond_ini = input('Quelle condition initiale ? 1) u(x,y) = 1.B(x,h) ');
    if cond_ini==1
        u0 = zeros(size(Centre_tetra,1),1);
        for i=1:size(Centre_tetra,1)
            if pdist([Centre_tetra(i,:);lambda*[0.25,0.25,0]])<5*h
                u0(i,1) = 1;
            end
        end
    end
    Delta_t =input('Quelle valeur pour Delta_t ? ');
end

%% Approximation terme de convection
if equation==2||equation==3
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
end

%% Création des fonctions f, g et h
if equation==1
    probleme = input('Quel exemple ? 1) u(x,y,z)=1 2) u(x,y,z) = x^2+y^2+z^2 3) u(x,y,z) = x+y+z ');
    if probleme==1
        f_dif = @(x) 0;
        g_dif = @(x) 1;
        h_dif = @(x) 0;
    elseif probleme==2
        f_dif = @(x) -6;
        g_dif = @(x) norm(x,2)^2;
        h_dif = @(x) 2*[x(1),x(2),x(3)]*[x(1);x(2);x(3)]/norm([x(1),x(2),x(3)],2);
    elseif probleme==3
        f_dif = @(x) 0;
        g_dif = @(x) sum(x);
        h_dif = @(x) x*[1;1;1]/norm(x,2);
    elseif probleme==4
        f_dif = @(x) 75*pi^2*cos(5*pi*sum(x));
        g_dif = @(x) cos(5*pi*sum(x));
        h_dif = @(x) -x*5*pi*sin(5*pi*sum(x))*[1;1;1]/norm(x,2);
    end
elseif equation==2
    f_conv = @(x) 0;
    g_conv = @(x) 0;
    h_conv = @(x) 0;
elseif equation ==3
    f_func = @(x) 0;
    g_func = @(x) 0;
    h_func = @(x) 0;
elseif equation==4
    nbr_en=10;
    r_alea = lambda*rand(nbr_en,1);
    theta_alea = pi*rand(nbr_en,1);
    phi_alea = 2*pi*rand(nbr_en,1);
    alea = r_alea.*[sin(theta_alea).*cos(phi_alea),sin(theta_alea).*sin(phi_alea),cos(theta_alea)];
    f_nut = @(x) sum(1*(pdist([x;alea])<5*h));
    g_nut = @(x) 0;
    h_nut = @(x) 0;
elseif equation==5
    nbr_en=10;
    r_alea = lambda*rand(nbr_en,1);
    theta_alea = pi*rand(nbr_en,1);
    phi_alea = 2*pi*rand(nbr_en,1);
    alea = r_alea.*[sin(theta_alea).*cos(phi_alea),sin(theta_alea).*sin(phi_alea),cos(theta_alea)];
    f_nut = @(x) sum(1*(pdist([x;alea])<5*h));
    g_nut = @(x) 0;
    h_nut = @(x) 0;
    f_cel = @(x) 0;
    g_cel = @(x) 0;
    h_cel = @(x) 0;
end
%% Création des matrices A et b et résolution de Au=b
if equation==1
    [A,b] = matrice_tetra(Tetra,f_dif,g_dif,h_dif,equation);
    u = A\b; % valeur au centre des triangles
    
elseif equation==2
    [A,b] = matrice_tetra(Tetra,f_conv,g_conv,h_conv,equation);
    T = input('T = ? ');
    i=1;
    figure('units','normalized','outerposition',[0 0 1 1]);
    u=u0;
    while (i*Delta_t<T)
        u = (speye(size(u,1))+Delta_t*A)*u +Delta_t*b;%Schéma explicite
        %u = (2*eye(size(u,1))-A)\(u+b);%Schéma implicite
        if i*Delta_t*10==floor(i*Delta_t*10)
            for k=1:3
                for l=1:5
                    subplot(3,5,(k-1)*5+l)
                    axis([-lambda lambda -lambda lambda]);
                    truc=(Centre_tetra(:,k)<lambda/4*(l-3)+l/2).*(Centre_tetra(:,k)>lambda/4*(l-3)-l/2);
                    B = Centre_tetra(truc==1,:);
                    if k==1
                        scatter(B(:,2),B(:,3),[],u(truc==1,1),'filled');
                        title(['Plan y-z avec x=',num2str(lambda/4*(l-3))]);
                    elseif k==2
                        scatter(B(:,1),B(:,3),[],u(truc==1,1),'filled');
                        title(['Plan x-z avec y=',num2str(lambda/4*(l-3))]);
                    else
                        scatter(B(:,1),B(:,2),[],u(truc==1,1),'filled');
                        title(['Plan x-y avec z=',num2str(lambda/4*(l-3))]);
                    end
                end
            end
            drawnow
        end
        i=i+1;
    end
    figure;
    scatter3(Centre_tetra(:,1),Centre_tetra(:,2),Centre_tetra(:,3),20,u,'filled')
    
elseif equation==3
    [A,b] = matrice_tetra(Tetra,f_func,g_func,h_func,equation);
    T = input('T = ? ');
    i=1;
    figure('units','normalized','outerposition',[0 0 1 1]);
    colormap(flipud(hot));
    colorbar;
    u=u0;
    while (i*Delta_t<T)
        u = (speye(size(u,1))+Delta_t*A)*u +Delta_t*b;%Schéma explicite
        %u = (2*eye(size(u,1))-A)\(u+b);%Schéma implicite
        if i*Delta_t*10==floor(i*Delta_t*10)
            colorbar;
            for k=1:3
                for l=1:5
                    subplot(3,5,(k-1)*5+l)
                    axis([-lambda lambda -lambda lambda]);
                    truc=(Centre_tetra(:,k)<lambda/4*(l-3)+l/2).*(Centre_tetra(:,k)>lambda/4*(l-3)-l/2);
                    B = Centre_tetra(truc==1,:);
                    if k==1
                        scatter(B(:,2),B(:,3),[],u(truc==1,1),'filled');
                        title(['Plan y-z avec x=',num2str(lambda/4*(l-3))]);
                    elseif k==2
                        scatter(B(:,1),B(:,3),[],u(truc==1,1),'filled');
                        title(['Plan x-z avec y=',num2str(lambda/4*(l-3))]);
                    else
                        scatter(B(:,1),B(:,2),[],u(truc==1,1),'filled');
                        title(['Plan x-y avec z=',num2str(lambda/4*(l-3))]);
                    end
                end
            end
            drawnow
        end
        i=i+1;
    end
    figure;
    scatter3(Centre_tetra(:,1),Centre_tetra(:,2),Centre_tetra(:,3),20,u,'filled')
    colormap(flipud(hot));
    colorbar;
    
elseif equation==4
    [A,b] = matrice_tetra(Tetra,f_nut,g_nut,h_nut,equation);
    u = A\b; % valeur au centre des triangles
    figure;
    scatter3(Centre_tetra(:,1),Centre_tetra(:,2),Centre_tetra(:,3),20,u,'filled')
    colormap(flipud(hot));
    axis square;
    colorbar;
    
elseif equation==5
    u=u0;
    gamma = 1;%input('Gamma = ');
    chi = 100;%input('chi = ');
    k=2;
    T = 100;%input('T = ');
    [A_nut,b_nut]=matrice_tetra(Tetra,f_nut,g_nut,h_nut,4);
    A_nut_act = A_nut;
    for i=1:size(A_nut,1)
        A_nut_act(i,i) = A_nut(i,i) + Volume(i,1)*gamma*u(i,1);
    end
    Nutriment = A_nut_act\b_nut;
    [A_cel,b_cel]=matrice_tetra(Tetra,f_cel,g_cel,h_cel,equation);    
    for i=1:size(A_cel,1)
        for j=1:size(TN,2)
            if isnan(TN(i,j))
            else
                A_cel(i,i) = A_cel(i,i) + chi*Coef_trans(i,j)*(Nutriment(TN(i,j))-Nutriment(i))*(Nutriment(TN(i,j))-Nutriment(i)>0);
                A_cel(i,TN(i,j)) = A_cel(i,TN(i,j)) + chi*Coef_trans(i,j)*(Nutriment(TN(i,j))-Nutriment(i))*(Nutriment(TN(i,j))-Nutriment(i)<0);
            end
        end
    end
    u = (speye(size(u,1))-Delta_t*A_cel)*u +Delta_t*b_cel;
    figure('units','normalized','outerposition',[0 0 1 1]);
    while (k*Delta_t<T)
        A_nut_act = A_nut;
        for i=1:size(A_nut,1)
            A_nut_act(i,i) = A_nut(i,i) + Volume(i,1)*gamma*u(i,1);
        end
        Nutriment = A_nut_act\b_nut;
        for i=1:size(A_cel,1)
            for j=1:size(TN,2)
                if isnan(TN(i,j))
                else
                A_cel(i,i) = A_cel(i,i) + chi*Coef_trans(i,j)*(Nutriment(TN(i,j))-Nutriment(i))*(Nutriment(TN(i,j))-Nutriment(i)>0);
                A_cel(i,TN(i,j)) = A_cel(i,TN(i,j)) + chi*Coef_trans(i,j)*(Nutriment(TN(i,j))-Nutriment(i))*(Nutriment(TN(i,j))-Nutriment(i)<0);
                end
            end
        end
        u = (speye(size(u,1))-Delta_t*A_cel)*u +Delta_t*b_cel;
        if k*Delta_t*50==floor(k*Delta_t*50)
            subplot(1,2,1);
            scatter3(Centre_tetra(:,1),Centre_tetra(:,2),Centre_tetra(:,3),20,u,'filled')
            axis square;
            alpha(.5)
            colormap(flipud(hot));
            colorbar;
            title('Tumor cells (concentration)');
            subplot(1,2,2);
            scatter3(Centre_tetra(:,1),Centre_tetra(:,2),Centre_tetra(:,3),20,Nutriment,'filled')
            axis square;
            alpha(.5)
            colorbar;
            title('Nutrients (concentration)');
            suptitle(['t = ',num2str(k*Delta_t),' s']);
            drawnow
        end
        k=k+1;
    end
    
end

%% Calcul des vraies valeurs
if equation==1
     z = zeros(size(Centre_tetra,1),1);
     for i=1:size(Centre_tetra,1)
         z(i) = g([Centre_tetra(i,1),Centre_tetra(i,2),Centre_tetra(i,3)]); 
     end
end

%% Les plots 
if equation==1
    figure;
    scatter3(Centre_tetra(:,1),Centre_tetra(:,2),Centre_tetra(:,3),10,u,'filled')
    figure;
    scatter3(Centre_tetra(:,1),Centre_tetra(:,2),Centre_tetra(:,3),10,z,'filled')
    disp(max(norm((u-z)./z,2)));
end