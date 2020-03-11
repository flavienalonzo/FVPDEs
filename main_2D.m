h=0.05;lambda=5;uniform=1;condition=2;
V = @(x) [-x(2),x(1)];
%% 
Tri = CreateTriMesh2D(h,lambda,uniform,condition,V);
save(['Tri_2D_',num2str(h),'_',num2str(lambda),'_',num2str(uniform),'_',num2str(condition),'_',func2str(V),'.mat'],'Tri');

%% Ecrire les équations 
disp('1) -D\nabla^{2}u(x) = f(x)');
disp('2) \partial_{t}u(t,x) + \vec V(x).\nabla u(t,x) = f(x)');
disp('3) \partial_{t}u(t,x) -D\nabla^{2}u(t,x) + \vec V(x).\nabla u(t,x) = f(x)');
disp('4) -D\nabla^{2}u(t,x) + \beta u(t,x) = f(x)');
disp('5) \partial_{t}u(t,x) -D\nabla^{2}u(t,x) + \beta u(t,x) = f(x)');

%% 
load(['Tri_2D_',num2str(h),'_',num2str(lambda),'_',num2str(uniform),'_',num2str(condition),'_',func2str(V),'.mat']);
equation=6;

%% 
Points = Tri{1};
ConnectivityList = Tri{2};
Eall = Tri{3};
Typesegment = Tri{4};
TN = Tri{5};
Centre_tri = Tri{6};
Dist_tri = Tri{7};
Coef_trans = Tri{8};
Volume = Tri{9};
Long_segment = Tri{10};

%% 
u0 = zeros(size(Centre_tri,2),1);
for i=1:size(Centre_tri,2)
    if pdist([Centre_tri(:,i)';lambda*[0.15,0.15]])<10*h
        u0(i,1) = 1;
    end
end
u_e0 = zeros(size(Centre_tri,2),1);
unif = linspace(0,2*pi,1000);
for i=1:size(Centre_tri,2)
    %if pdist([Centre_tri(:,i)';lambda*[-0.25,0.25]])<100*h
    lissajous = pdist([Centre_tri(:,i)';sqrt(lambda)*[sin(5*unif'),cos(3*unif')]]);
    u_e0(i,1) = sum(lissajous(1,1:size(unif,2))<h);
        %= exp(-pdist([Centre_tri(:,i)';lambda*[-0.25,0.25]]))+exp(-pdist([Centre_tri(:,i)';lambda*[0.25,-0.25]]));
    %end
end
Delta_t = 0.01;
%% 
format long;
h=0.01;
fd = @(p) sqrt(sum(p.^2,2))-1;
fh = @(p) max(5*sqrt(sum(p.^2,2))-1,2);%ones(size(p,1),1);
[p,t] = distmesh(fd,fh,h,[-1,-1;1,1],[],'IT_MAX',100000);%@huniform
% save('maillage.mat','p');

%% 
c = 5;
dt = delaunayTriangulation(c*p);
ch = convexHull(dt);
condition = input('Condition au bord ? 1 - Dirichlet, 2 - Neumann ');
type(ch)= condition;
x = p(:,1);y=p(:,2);
% figure;
% triplot(dt);
NbTri = size(dt.ConnectivityList,1);

%% Création de la liste des segments du maillage
[E,Eall]=Edge(dt);
NbSeg = size(Eall,1);

%% Attribution des segments à chaque triangle
Tri={};
Tri{1} = dt.Points;
Tri{2} = dt.ConnectivityList;
Tri{3} = Eall;

%% Ajout de la longueur des segments
Long_segment = zeros(1,NbSeg);
for i=1:NbSeg
    Long_segment(1,i) = pdist(dt.Points(Eall(i,:),:),'euclidean');
end
Seg{3} = Long_segment;
Tri{10} = Long_segment;

%% Check if all triangles are acute
NewPoints = [];
for i=1:NbTri
    if sum(Long_segment(1,3*(i-1)+1:3*(i-1)+3).^2)<=2*max(Long_segment(1,3*(i-1)+1:3*(i-1)+3))^2
        if (Long_segment(1,3*(i-1)+1)==max(Long_segment(1,3*(i-1)+1:3*(i-1)+3)))
            k=1;
            NewPoints(end+1,:) = mean(dt.Points(dt.ConnectivityList(i,[1,2]),:))./c;
        elseif (Long_segment(1,3*(i-1)+2)==max(Long_segment(1,3*(i-1)+1:3*(i-1)+3)))
            k=2;
            NewPoints(end+1,:) = mean(dt.Points(dt.ConnectivityList(i,[2,3]),:))./c;
        else 
            k=3;
            NewPoints(end+1,:) = mean(dt.Points(dt.ConnectivityList(i,[1,3]),:))./c;
        end
        
    end
end

%% La Liste des triangles voisins d'un triangle
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

%% Ajout du type de segment
%0 - segment à l'interieur du domaine
%1 - segment de type Dirichlet
%2 - segment de type Neumann
%ch = convexHull(dt);
Seg = {};
Seg{1} = Eall;
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
Seg{2} = Typesegment; Tri{4} = Typesegment; 


%% Ajout des cordonnées de chaque triangle
Tri{5} = TN;
Centre_tri=circumcenter(dt)';
Tri{6} = Centre_tri;
% figure;
% plot(Centre_tri(1,:),Centre_tri(2,:),'+')

%% Ajout des distances entre triangles
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

%% Création du coefficient de transmissibilité
Coef_trans = 0*Dist_tri;
for i=1:NbTri
    Coef_trans(i,1) = Long_segment(1,3*i-2)/Dist_tri(i,1) ;
    Coef_trans(i,2) = Long_segment(1,3*i-1)/Dist_tri(i,2) ;
    Coef_trans(i,3) = Long_segment(1,3*i)/Dist_tri(i,3) ;
end

Tri{8} = Coef_trans;

%% Ajout du volume de chaque triangle
Volume = zeros(1,NbTri);
for i=1:NbTri
    s = 0.5*sum(Long_segment(3*i-2:3*i));
    Volume(1,i) = sqrt(s*(s-Long_segment(3*i-2))*(s-Long_segment(3*i-1))*(s-Long_segment(3*i)));
end
Tri{9} = Volume;

%% Introduction au système à résoudre
equation = input('Quelle équation ? 1) Diffusion 2) Convection 3) Diffusion-convection ');
Delta_t=0;
if equation==2 || equation==3
    vitesse = input('Quelle vitesse de convection ? 1) V=(-y,x) 2) V=(-y,x)/norm(x,y) 3) V=(-y,x)/norm(1.5*x,y) ');
    if vitesse==1
        V = @(x) [-x(2),x(1)];
    elseif vitesse==2
        V = @(x) [-x(2),x(1)]./norm(x,2);
    elseif vitesse==3
        V = @(x) [-x(2),x(1)]./(norm(x,2)+0.5*x(1)^2);
    end
    cond_ini = input('Quelle condition initiale ? 1) u(x,y) = 1.B(x,h) ');
    if cond_ini==1
        u0 = zeros(size(Centre_tri,2),1);
        for i=1:size(Centre_tri,2)
            if pdist([Centre_tri(:,i)';lambda*[0.25,0.25]])<2*h
                u0(i,1) = 1;
            end
        end
    end
    Delta_t =input('Quelle valeur pour Delta_t ? ');
end

%% Approximation d'un terme de convection
if equation==2 || equation==3
    Vites_seg = zeros(NbSeg,2);
    Normales_seg = zeros(NbSeg,2);
    Prodvitnor_seg = zeros(NbSeg,1);
    for i=1:NbSeg
        Vites_seg(i,:) = V(0.5*(dt.Points(Eall(i,1),:)+dt.Points(Eall(i,2),:)));
        %Vites_seg(i,1)=-0.5*(dt.Points(Eall(i,1),2)+dt.Points(Eall(i,2),2));%-1/Long_segment(1,i)*(dt.Points(Eall(i,2),1)-dt.Points(Eall(i,1),1))*(dt.Points(Eall(i,2),2)^2-dt.Points(Eall(i,1),2)^2);
        %Vites_seg(i,2)=0.5*(dt.Points(Eall(i,1),1)+dt.Points(Eall(i,2),1));%1/Long_segment(1,i)*(dt.Points(Eall(i,2),2)-dt.Points(Eall(i,1),2))*(dt.Points(Eall(i,2),1)^2-dt.Points(Eall(i,1),1)^2);
%             ax = (dt.Points(Eall(i,2),2)-dt.Points(Eall(i,1),2))/(dt.Points(Eall(i,2),1)-dt.Points(Eall(i,1),1));
%             Normales_seg(i,2) = ax/norm([ax,1],2);
%             Normales_seg(i,1) = 1/norm([ax,1],2);
        Normales_seg(i,:) = (0.5*(dt.Points(Eall(i,1),:)+dt.Points(Eall(i,2),:))-Centre_tri(:,floor((i-1)/3)+1)')/(norm((0.5*(dt.Points(Eall(i,1),:)+dt.Points(Eall(i,2),:))-Centre_tri(:,floor((i-1)/3)+1)'),2));
%             if Normales_seg(i,:)'*GM<0
%                 Normales_seg(i,:) = -Normales_seg(i,:);
%             end
        Prodvitnor_seg(i,1) = -Vites_seg(i,:)*Normales_seg(i,:)';
    end
    Tri{11} = Vites_seg;
    Tri{12} = Normales_seg;
    Tri{13} = Prodvitnor_seg;
end

%% Création des fonctions f, g et h
if equation==1
    probleme = input('Quel probleme ? 1) Quadratique 2) u(x,y)=1 3) u(x,y)=x+y 4) u(x,y)=x^2-y^2 5) u(x,y)=cos(5pi(x+y)) 6) u(x,y)=0 ');
    if probleme==1
        f_dif = @(x) -2;
        g_dif = @(x) 0.5*x*[1,0.5;0.5,1]*x'-x*[1;-1]+0.6;
        h_dif = @(x) x*([1,0.5;0.5,1]*x'-[1;-1])/norm(x,2);
    elseif probleme==2
        f_dif = @(x) 0;
        g_dif = @(x) 1;
        h_dif = @(x) 0;
    elseif probleme==3
        f_dif = @(x) 0;
        g_dif = @(x) sum(x);
        h_dif = @(x) x*[1;1]/norm(x,2);
    elseif probleme==4
        f_dif = @(x) 0;
        g_dif = @(x) x(1)^2-x(2)^2;
        h_dif = @(x) x*[2*x(1);-2*x(2)]/norm(x,2);
    elseif probleme==5
        f_dif = @(x) 2*25*pi^2*cos(5*pi*sum(x));
        g_dif = @(x) cos(5*pi*sum(x));
        h_dif = @(x) -x*5*pi*sin(5*pi*sum(x))*[1;1]/norm(x,2);
    elseif probleme==6
        f_dif = @(x) 0;
        g_dif = @(x) 0;
        h_dif = @(x) 0;
    end
end
if equation==2
    f_conv = @(x) 0;
    g_conv = @(x) 0;
    h_conv = @(x) 0;
elseif equation==3
    f = @(x) 0;
    g = @(x) 0;
    h = @(x) 0;
elseif equation==4
    nbr_en=10;
    r_alea = lambda*rand(nbr_en,1);
    theta_alea = 2*pi*rand(nbr_en,1);
    alea = r_alea.*[cos(theta_alea),sin(theta_alea)];
    f_nut = @(x) sum(1*(pdist([x;alea])<5*h));
    g_nut = @(x) 0;
    h_nut = @(x) 0;
elseif equation==5
%     nbr_en=50;
%     r_alea = lambda*rand(nbr_en,1);
%     theta_alea = 2*pi*rand(nbr_en,1);
%     alea = r_alea.*[cos(theta_alea),sin(theta_alea)];
    f_nut = @(x) 0;%sum(1*(pdist([x;alea])<3*h));
    g_nut = @(x) 0;
    h_nut = @(x) 0;
    f_cel = @(x) 0;
    g_cel = @(x) 0;
    h_cel = @(x) 0;
elseif equation==6
    f_nut = @(x) 0;
    g_nut = @(x) 0;
    h_nut = @(x) 0;
    f_cel = @(x) 0;
    g_cel = @(x) 0;
    h_cel = @(x) 0;
    f_end = @(x) 0;
    g_end = @(x) 0;
    h_end = @(x) 0;
end

%% Création des matrices A et b et résolution de Au=b
if equation==1
    [A,b] = matrice_tri(Tri,f_dif,g_dif,h_dif,equation);
    u = A\b; % valeur au centre des triangles
    
elseif equation==2
    [A,b] = matrice_tri(Tri,f_conv,g_conv,h_conv,equation);
    schema = input('Schéma : 1) explicite 2) implicite ? ');
    Delta_t = input('\Delta t = ');
    T = input('T = ? ');
    i=1;
    figure;
    u=u0;
    scatter(Centre_tri(1,:),Centre_tri(2,:),[],u0,'filled')
    title(['t = ',num2str(0),' s']);
    axis square;
    while (i*Delta_t<T)
        if schema==1
            u = (speye(size(u,1))-Delta_t*A_cel./Volume)*u +Delta_t*b_cel./Volume';
        elseif schema==2
            u = (speye(size(u,1))+Delta_t*A_cel./Volume)\(u +Delta_t*b_cel./Volume');
        end
        if i*Delta_t*10==floor(i*Delta_t*10)
            scatter(Centre_tri(1,:),Centre_tri(2,:),[],u,'filled')
            axis square;
            title(['t = ',num2str(i*Delta_t),' s']);
            drawnow
        end
        i=i+1;
    end
    scatter(Centre_tri(1,:),Centre_tri(2,:),[],u,'filled')
    axis square;
    title(['t = ',num2str(i*Delta_t),' s']);
    drawnow
    
elseif equation==3
    [A,b] = matrice_tri(Tri,f,g,h,equation);
    schema = input('Schéma : 1) explicite 2) implicite ? ');
    Delta_t = input('\Delta t = ');
    T = input('T = ? ');
    i=1;
    figure('units','normalized','outerposition',[0 0 1 1]);
    axis square;
    colormap(flipud(hot));
    colorbar;
    u=u0;
    scatter(Centre_tri(1,:),Centre_tri(2,:),[],u0,'filled')
    title(['t = ',num2str(0),' s']);
    while (i*Delta_t<T)
        if schema==1
            u = (speye(size(u,1))-Delta_t*A_cel./Volume)*u +Delta_t*b_cel./Volume';
        elseif schema==2
            u = (speye(size(u,1))+Delta_t*A_cel./Volume)\(u +Delta_t*b_cel./Volume');
        end
        if i*Delta_t*10==floor(i*Delta_t*10)
            scatter(Centre_tri(1,:),Centre_tri(2,:),[],u,'filled')
            axis square;
            colorbar;
            title(['t = ',num2str(i*Delta_t),' s']);
            drawnow
        end
        i=i+1;
    end
    scatter(Centre_tri(1,:),Centre_tri(2,:),[],u,'filled')
    axis square;
    colorbar;
    title(['t = ',num2str(i*Delta_t),' s']);
    drawnow
    
elseif equation==4
    [A,b] = matrice_tri(Tri,f_nut,g_nut,h_nut,equation);
    u = A\b; % valeur au centre des triangles
    figure;
    scatter(Centre_tri(1,:),Centre_tri(2,:),5,u,'filled')
    colormap(flipud(hot));
    axis square;
    colorbar;
    
elseif equation==5
    u=u0;
    gamma = 1;%input('Gamma = ');
    chi = 100;%input('chi = ');
    schema = input('Schéma : 1) explicite 2) implicite ? ');
    alpha = 1;
    k=2;
    tumorsize = [Volume*(u0>10e-5)];
    tumormass = [Volume*u0];
    T = 100;%input('T = ');
    [A_nut,b_nut]=matrice_tri(Tri,f_nut,g_nut,h_nut,4);
    A_nut_act = A_nut;
    for i=1:size(A_nut,1)
        A_nut_act(i,i) = A_nut(i,i) + Volume(1,i)*gamma*u(i,1);
        b_nut_act(i,1) = b_nut(i,1) + Volume(1,i)*alpha*u_e0(i,1);
    end
    Nutriment = A_nut_act\b_nut_act;
    [A_cel,b_cel]=matrice_tri(Tri,f_cel,g_cel,h_cel,equation);    
    for i=1:size(A_cel,1)
        for j=1:size(TN,2)
            if isnan(TN(i,j))
            else
                A_cel(i,i) = A_cel(i,i) + chi*Coef_trans(i,j)*(Nutriment(TN(i,j))-Nutriment(i))*(Nutriment(TN(i,j))-Nutriment(i)>0);
                A_cel(i,TN(i,j)) = A_cel(i,TN(i,j)) + chi*Coef_trans(i,j)*(Nutriment(TN(i,j))-Nutriment(i))*(Nutriment(TN(i,j))-Nutriment(i)<0);
            end
        end
    end
    if schema==1
        Delta_t = input('\Delta t = ');
        u = (speye(size(u,1))-Delta_t*A_cel./Volume)*u +Delta_t*b_cel./Volume';
    elseif schema==2
        Delta_t = input('\Delta t = ');
        u = mldivide(speye(size(u,1))+Delta_t*A_cel./Volume,u +Delta_t*b_cel./Volume');
    end
    tumorsize(end+1) = Volume*(u>10e-5);
    tumormass(end+1) = Volume*u;
    figure('units','normalized','outerposition',[0 0 1 1]);
    while (k*Delta_t<T)
        A_nut_act = A_nut;
        for i=1:size(A_nut,1)
            A_nut_act(i,i) = A_nut(i,i) + Volume(1,i)*gamma*u(i,1);
            b_nut_act(i,1) = b_nut(i,1) + Volume(1,i)*alpha*u_e0(i,1);
        end
        Nutriment = A_nut_act\b_nut_act;
        for i=1:size(A_cel,1)
            for j=1:size(TN,2)
                if isnan(TN(i,j))
                else
                A_cel(i,i) = A_cel(i,i) + chi*Coef_trans(i,j)*(Nutriment(TN(i,j))-Nutriment(i))*(Nutriment(TN(i,j))-Nutriment(i)>0);
                A_cel(i,TN(i,j)) = A_cel(i,TN(i,j)) + chi*Coef_trans(i,j)*(Nutriment(TN(i,j))-Nutriment(i))*(Nutriment(TN(i,j))-Nutriment(i)<0);
                end
            end
        end
        if schema==1
            u = (speye(size(u,1))-Delta_t*A_cel./Volume)*u +Delta_t*b_cel./Volume';
        elseif schema==2
            u = mldivide(speye(size(u,1))+Delta_t*A_cel./Volume,u +Delta_t*b_cel./Volume');
        end
        tumorsize(end+1) = Volume*(u>10e-5);
        tumormass(end+1) = Volume*u;
        if k*Delta_t*10==floor(k*Delta_t*10)
            ax(1)=subplot(2,2,1);
            scatter(Centre_tri(1,:),Centre_tri(2,:),2.5,u,'filled')
            colormap(ax(1),flipud(hot));
            axis square;
            c1 = colorbar;
            c1.Label.String = 'u (concentration)';
            title({'Tumor cells','\partial_{t}u-D_{u}\nabla^{2}u +\chi\nabla.(u\nabla c)= ru '});
            ax(2)=subplot(2,2,2);
            colormap(ax(2),flipud(summer));
            scatter(Centre_tri(1,:),Centre_tri(2,:),2.5,Nutriment,'filled')
            axis square;
            c2 = colorbar;
            c2.Label.String = 'c (concentration)';
            title({'Nutrients','-D_{c}\nabla^{2} c+\beta c = \alpha u_{e} - \gamma u c'});
            ax(3)=subplot(2,2,[3,4]);
            hold on;
            yyaxis left;
            plot(0:Delta_t:k*Delta_t,tumorsize);
            yyaxis right;
            plot(0:Delta_t:k*Delta_t,tumormass);
            legend({'Tumor size','Tumor mass'},'Location','northwest');
            xlabel('Time');
            hold off;
            title({'Tumor volume and mass'});
            %
            suptitle({['t = ',num2str(k*Delta_t)],['\Delta t = ',num2str(Delta_t)]});
            drawnow
        end
        k=k+1;
    end
    
    elseif equation==6
    u=u0;
    u_e=u_e0;
    alpha = 1;%input('Alpha = ');
    gamma = 1;%input('Gamma = ');
    chi = 1;%input('chi = ');
    schema = input('Schéma : 1) explicite 2) implicite ? ');
    tumorsize = [Volume*(u0>10e-5)];
    tumormass = [Volume*u0];
    k=2;
    T = 100;%input('T = ');
    [A_nut,b_nut]=matrice_tri(Tri,f_nut,g_nut,h_nut,4);
    A_nut_act = A_nut;
    b_nut_act = b_nut;
    for i=1:size(A_nut,1)
        %if alpha*u_e(i,1)>gamma*u(i,1)
            A_nut_act(i,i) = A_nut(i,i) + Volume(1,i)*gamma*u(i,1);
            b_nut_act(i,1) = b_nut(i,1) + Volume(1,i)*alpha*u_e(i,1);
        %end
    end
    Nutriment = A_nut_act\b_nut_act;
    [A_cel,b_cel]=matrice_tri(Tri,f_cel,g_cel,h_cel,5);    
    for i=1:size(A_cel,1)
        for j=1:size(TN,2)
            if isnan(TN(i,j))
            else
                A_cel(i,i) = A_cel(i,i) + chi*Coef_trans(i,j)*(Nutriment(TN(i,j))-Nutriment(i))*(Nutriment(TN(i,j))-Nutriment(i)>0);
                A_cel(i,TN(i,j)) = A_cel(i,TN(i,j)) + chi*Coef_trans(i,j)*(Nutriment(TN(i,j))-Nutriment(i))*(Nutriment(TN(i,j))-Nutriment(i)<0);
            end
        end
    end
    if schema==1
        Delta_t = input('\Delta t = ');
        u = (speye(size(u,1))-Delta_t*A_cel./Volume)*u +Delta_t*b_cel./Volume';
    elseif schema==2
        Delta_t = input('\Delta t = ');
        u = (speye(size(u,1))+Delta_t*A_cel./Volume)\(u +Delta_t*b_cel./Volume');
    end
    [A_end,b_end]=matrice_tri(Tri,f_end,g_end,h_end,5);
    if schema==1
        u_e = (speye(size(u_e,1))-Delta_t*A_end./Volume)*u_e +Delta_t*b_end./Volume';
    elseif schema==2
        u_e = (speye(size(u_e,1))+Delta_t*A_end./Volume)\(u_e +Delta_t*b_end./Volume');
    end
    tumorsize(end+1) = Volume*(u>10e-5);
    tumormass(end+1) = Volume*u;    
    figure('units','normalized','outerposition',[0 0 1 1]);
    while (k*Delta_t<T)
        A_nut_act = A_nut;
        for i=1:size(A_nut,1)
           % if alpha*u_e(i,1)>gamma*u(i,1)
                A_nut_act(i,i) = A_nut(i,i) + Volume(1,i)*gamma*u(i,1);
                b_nut_act(i,1) = b_nut(i,1) + Volume(1,i)*alpha*u_e(i,1);
            %end
        end
        Nutriment = A_nut_act\b_nut_act;
        for i=1:size(A_cel,1)
            for j=1:size(TN,2)
                if isnan(TN(i,j))
                else
                A_cel(i,i) = A_cel(i,i) + chi*Coef_trans(i,j)*(Nutriment(TN(i,j))-Nutriment(i))*(Nutriment(TN(i,j))-Nutriment(i)>0);
                A_cel(i,TN(i,j)) = A_cel(i,TN(i,j)) + chi*Coef_trans(i,j)*(Nutriment(TN(i,j))-Nutriment(i))*(Nutriment(TN(i,j))-Nutriment(i)<0);
                end
            end
        end
        if schema==1
            u = (speye(size(u,1))-Delta_t*A_cel./Volume)*u +Delta_t*b_cel./Volume';
        elseif schema==2
            u = (speye(size(u,1))+Delta_t*A_cel./Volume)\(u +Delta_t*b_cel./Volume');
        end
        if schema==1
            u_e = (speye(size(u_e,1))-Delta_t*A_end./Volume)*u_e +Delta_t*b_end./Volume';
        elseif schema==2
            u_e = (speye(size(u_e,1))+Delta_t*A_end./Volume)\(u_e +Delta_t*b_end./Volume');
        end
        tumorsize(end+1) = Volume*(u>10e-5);
        tumormass(end+1) = Volume*u;
        if k*Delta_t*50000==floor(k*Delta_t*50000)
            ax(1)=subplot(2,3,1);
            scatter(Centre_tri(1,:),Centre_tri(2,:),2.5,u,'filled')
            colormap(ax(1),flipud(hot));
            axis square;
            grid(ax(1),'on');
            grid minor;
            c1 = colorbar;
            c1.Label.String = 'u (concentration)';
            title({'Tumor cells','\partial_{t}u-D_{u}\nabla^{2}u +\chi\nabla.(u\nabla c)= ru '});
            ax(2)=subplot(2,3,2);
            colormap(ax(2),flipud(summer));
            scatter(Centre_tri(1,:),Centre_tri(2,:),2.5,Nutriment,'filled')
            axis square;
            grid(ax(2),'on');
            grid minor;
            c2 = colorbar;
            c2.Label.String = 'c (concentration)';
            title({'Nutrients',' -D_{c}\nabla^{2} c+\beta c = \alpha u_{e} - \gamma u c'});
            ax(3)=subplot(2,3,3);
            colormap(ax(3),flipud(bone));
            scatter(Centre_tri(1,:),Centre_tri(2,:),2.5,u_e,'filled')
            axis square;
            grid(ax(3),'on');
            grid minor;
            c3 = colorbar;
            c3.Label.String = 'u_{e} (concentration)';
            title({'Endothelial cells','\partial_{t}u_{e}-D_{e}\nabla^{2}u_{e} = r_{e}u_{e} '});
            ax(4)=subplot(2,3,[4,5,6]);
            hold on;
            yyaxis left;
            plot(0:Delta_t:k*Delta_t,tumorsize);
            yyaxis right;
            plot(0:Delta_t:k*Delta_t,tumormass);
            legend({'Tumor size','Tumor mass'},'Location','northwest');
            xlabel('Time');
            hold off;
            title({'Tumor volume and mass'});
            %
            suptitle({['t = ',num2str(k*Delta_t),' , ','\Delta t = ',num2str(Delta_t)]});
            drawnow
        end
        k=k+1;
    end
end

%% Calcul des vraies valeurs
if equation==1
     z = zeros(size(x));
     for i=1:NbTri
         z(i) = g([Centre_tri(1,i),Centre_tri(2,i)]); 
     end
end

%% Graphiques
if equation==1
    figure;
    hold on;
    %triplot(dt);
    plot3(Centre_tri(1,:),Centre_tri(2,:),u,'+')
    plot3(Centre_tri(1,:),Centre_tri(2,:),z,'.')
    hold off;
    disp(max(norm((u-z)./z,2)));
end

%% 
%figure;
%surf(X,Y,reshape(quad([x,y]),size(X,1),size(X,2)));

