function [A,b] = matrice_tetra(Tetra,f,g,h,equation)
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
    NbTetra = size(ConnectivityList,1);
    format long; 
    
    if equation==1
        A = sparse(NbTetra,NbTetra);
        b = sparse(NbTetra,1);
        for i=1:NbTetra       
           A(i,i) = sum(Coef_trans(i,:));
           b(i,1) = Volume(i,1)*f(Centre_tetra(i,:));
           for j=1:size(TN,2)
               if isnan(TN(i,j))
                   if Typetriang(4*(i-1)+j,1)==0

                   elseif Typetriang(4*(i-1)+j,1)==1
                       b(i,1) = b(i,1) + Coef_trans(i,j)*g(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   elseif Typetriang(4*(i-1)+j,1)==2
                       A(i,i) = A(i,i) - Coef_trans(i,j)+ Surftriang(4*(i-1)+j);
                       b(i,1) = b(i,1) + Surftriang(4*(i-1)+j)*h(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   end

               else
                   A(i,TN(i,j)) = -Coef_trans(i,j);
                   if Typetriang(4*(i-1)+j,1)==0

                   elseif Typetriang(4*(i-1)+j,1)==1
                       b(i,1) = b(i,1) + Coef_trans(i,j)*g(Points(TN(i,j),:));
                   elseif Typetriang(4*(i-1)+j,1)==2
                       A(i,i) = A(i,i) - Coef_trans(i,j)+ Surftriang(4*(i-1)+j);
                       b(i,1) = b(i,1) + Surftriang(4*(i-1)+j)*h(Points(TN(i,j),:));
                   end
               end
           end
        end
    elseif equation==2
        A = sparse(NbTetra,NbTetra);
        b = sparse(NbTetra,1);
        Prodvitnor_surf = Tetra{13};
        for i=1:NbTetra       
           b(i,1) = Volume(i,1)*f(Centre_tetra(i,:));
           for j=1:size(TN,2)
               if isnan(TN(i,j))
                   if Typetriang(4*(i-1)+j,1)==0
                       A(i,i) = A(i,i) - 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0);
                   elseif Typetriang(4*(i-1)+j,1)==1
                       A(i,i) = A(i,i) - 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0);
                       b(i,1) = b(i,1) + 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0)*g(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   end

               else
                   A(i,TN(i,j)) = 1/Volume(i,1)*Surftriang(4*(i-1)+j,1)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0);
                   if Typetriang(4*(i-1)+j,1)==0
                       A(i,i) = A(i,i) - 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0);
                   elseif Typetriang(4*(i-1)+j,1)==1
                       A(i,i) = A(i,i) - 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0);
                       b(i,1) = b(i,1) + 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0)*g(Points(TN(i,j),:));
                   end
               end
           end
        end
    elseif equation==3
        A = sparse(NbTetra,NbTetra);
        b = sparse(NbTetra,1);
        Prodvitnor_surf = Tetra{13};
        for i=1:NbTetra       
           A(i,i) = 1/1000000*sum(Coef_trans(i,:));
           b(i,1) = Volume(i,1)*f(Centre_tetra(i,:));
           for j=1:size(TN,2)
               if isnan(TN(i,j))
                   if Typetriang(4*(i-1)+j,1)==0
                       A(i,i) = A(i,i) - 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0);
                   elseif Typetriang(4*(i-1)+j,1)==1
                       A(i,i) = A(i,i) - 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0);
                       b(i,1) = b(i,1) + 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0)*g(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                       b(i,1) = b(i,1) + 1/1000000*Coef_trans(i,j)*g(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   elseif Typetriang(4*(i-1)+j,1)==2
                       A(i,i) = A(i,i) - 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0);
                       A(i,i) = A(i,i) - 1/1000000*Coef_trans(i,j)+ Surftriang(4*(i-1)+j);
                       b(i,1) = b(i,1) + Surftriang(4*(i-1)+j)*h(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   end

               else
                   A(i,TN(i,j)) = -1/1000000*Coef_trans(i,j) + 1/Volume(i,1)*Surftriang(4*(i-1)+j,1)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0);
                   if Typetriang(4*(i-1)+j,1)==0
                       A(i,i) = A(i,i) - 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0);
                   elseif Typetriang(4*(i-1)+j,1)==1
                       A(i,i) = A(i,i) - 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0);
                       b(i,1) = b(i,1) + 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0)*g(Points(TN(i,j),:));
                       b(i,1) = b(i,1) + 1/1000000*Coef_trans(i,j)*g(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   elseif Typetriang(4*(i-1)+j,1)==2
                       A(i,i) = A(i,i) - 1/Volume(i,1)*Surftriang(4*(i-1)+j)*Prodvitnor_surf(4*(i-1)+j,1)*(Prodvitnor_surf(4*(i-1)+j,1)>0);
                       A(i,i) = A(i,i) - 1/1000000*Coef_trans(i,j)+ Surftriang(4*(i-1)+j);
                       b(i,1) = b(i,1) + Surftriang(4*(i-1)+j)*h(Points(TN(i,j),:));
                   end
               end
           end
        end
    elseif equation==4
        beta = input('beta = ');
        D = input('D = ');
        A = sparse(NbTetra,NbTetra);
        b = sparse(NbTetra,1);
        wbm = waitbar(0,'Building A and b');
        for i=1:NbTetra       
            waitbar(i/NbTetra,wbm);
           A(i,i) = D*sum(Coef_trans(i,:))+beta*Volume(i,1);
           b(i,1) = Volume(i,1)*f(Centre_tetra(i,:));
           for j=1:size(TN,2)
               if isnan(TN(i,j))
                   if Typetriang(4*(i-1)+j,1)==0

                   elseif Typetriang(4*(i-1)+j,1)==1
                       b(i,1) = b(i,1) + D*Coef_trans(i,j)*g(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   elseif Typetriang(4*(i-1)+j,1)==2
                       A(i,i) = A(i,i) - D*Coef_trans(i,j);
                       b(i,1) = b(i,1) + D*Surftriang(4*(i-1)+j)*h(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');%'*([0 -1;1 0]*(Points(Eall(3*(i-1)+j,2),:)-Points(Eall(3*(i-1)+j,1),:))')/norm(Points(Eall(3*(i-1)+j,2),:)-Points(Eall(3*(i-1)+j,1),:),2);
                   end

               else
                   A(i,TN(i,j)) = -D*Coef_trans(i,j);
                   if Typetriang(4*(i-1)+j,1)==0

                   elseif Typetriang(4*(i-1)+j,1)==1
                       b(i,1) = b(i,1) + D*Coef_trans(i,j)*g(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   elseif Typetriang(4*(i-1)+j,1)==2
                       A(i,i) = A(i,i) - D*Coef_trans(i,j) ;
                       b(i,1) = b(i,1) + D*Surftriang(4*(i-1)+j)*h(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   end
               end
           end
        end
    elseif equation==5
        A = sparse(NbTetra,NbTetra);
        b = sparse(NbTetra,1);
        D = 1.05e-3;
        r = 1.07e-2;
        wbm = waitbar(0,'Building A and b');
        for i=1:NbTetra       
           waitbar(i/NbTetra,wbm);
           A(i,i) = D*sum(Coef_trans(i,:))-r*Volume(i,1);
           b(i,1) = Volume(i,1)*f(Centre_tetra(i,:));
           for j=1:size(TN,2)
               if isnan(TN(i,j))
                   if Typetriang(4*(i-1)+j,1)==0

                   elseif Typetriang(4*(i-1)+j,1)==1
                       b(i,1) = b(i,1) + D*Coef_trans(i,j)*g(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   elseif Typetriang(4*(i-1)+j,1)==2
                       A(i,i) = A(i,i) - D*Coef_trans(i,j);
                       b(i,1) = b(i,1) + D*Surftriang(4*(i-1)+j)*h(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   end

               else
                   A(i,TN(i,j)) = -D*Coef_trans(i,j);
                   if Typetriang(4*(i-1)+j,1)==0

                   elseif Typetriang(4*(i-1)+j,1)==1
                       b(i,1) = b(i,1) + D*Coef_trans(i,j)*g(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   elseif Typetriang(4*(i-1)+j,1)==2
                       A(i,i) = A(i,i) - D*Coef_trans(i,j) ;
                       b(i,1) = b(i,1) + D*Surftriang(4*(i-1)+j)*h(center(Points(Triang(4*(i-1)+j,1),:),Points(Triang(4*(i-1)+j,2),:),Points(Triang(4*(i-1)+j,3),:),'circumcenter')');
                   end
               end
           end
        end
    end
    close(wbm);
end