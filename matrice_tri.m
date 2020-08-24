function [A,b] = matrice_tri(Tri,f,g,h,equation)
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

    NbTri = size(ConnectivityList,1);
    format long; 
    if equation==1 
        A = sparse(NbTri,NbTri);
        b = sparse(NbTri,1);
        wbm = waitbar(0,'Building A and b');
        for i=1:NbTri       
            waitbar(i/NbTri,wbm);
           A(i,i) = sum(Coef_trans(i,:));
           b(i,1) = Volume(1,i)*f(Centre_tri(:,i)');
           for j=1:size(TN,2)
               if isnan(TN(i,j))
                   if Typesegment(1,3*(i-1)+j)==0

                   elseif Typesegment(1,3*(i-1)+j)==2
                       b(i,1) = b(i,1) + Coef_trans(i,j)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   elseif Typesegment(1,3*(i-1)+j)==1
                       A(i,i) = A(i,i) - Coef_trans(i,j);
                       b(i,1) = b(i,1) + Long_segment(1,3*(i-1)+j)*h(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));%'*([0 -1;1 0]*(Points(Eall(3*(i-1)+j,2),:)-Points(Eall(3*(i-1)+j,1),:))')/norm(Points(Eall(3*(i-1)+j,2),:)-Points(Eall(3*(i-1)+j,1),:),2);
                   end

               else
                   A(i,TN(i,j)) = -Coef_trans(i,j);
                   if Typesegment(1,3*(i-1)+j)==0

                   elseif Typesegment(1,3*(i-1)+j)==2
                       b(i,1) = b(i,1) + Coef_trans(i,j)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   elseif Typesegment(1,3*(i-1)+j)==1
                       A(i,i) = A(i,i) - Coef_trans(i,j) ;
                       b(i,1) = b(i,1) + Long_segment(1,3*(i-1)+j)*h(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   end
               end
           end
        end
    elseif equation==2
        Vites_seg = Tri{11};
        Normales = Tri{12};
        Prodvitnor_seg = Tri{13};
        A = sparse(NbTri,NbTri);
        b = sparse(NbTri,1);
        wbm = waitbar(0,'Building A and b');
        for i=1:NbTri       
            waitbar(i/NbTri,wbm);
           b(i,1) = Volume(1,i)*f(Centre_tri(:,i)');
           for j=1:size(TN,2)
               if isnan(TN(i,j))
                   if Typesegment(1,3*(i-1)+j)==0
                       A(i,i) = A(i,i) + Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1) - 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0);
                   else 
                       A(i,i) = A(i,i) + Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1) - 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0); 
                       b(i,1) = b(i,1) + 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   end

               else
                   A(i,TN(i,j)) = 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0);
                   if Typesegment(1,3*(i-1)+j)==0
                       A(i,i) = A(i,i) + Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1) - 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0);
                   else
                       A(i,i) = A(i,i) + Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1) - 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0);
                       b(i,1) = b(i,1) + Coef_trans(i,j)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:))) + 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   end
               end
           end
        end
    elseif equation==3
        A = sparse(NbTri,NbTri);
        b = sparse(NbTri,1);
        Prodvitnor_seg = Tri{13};
        wbm = waitbar(0,'Building A and b');
        for i=1:NbTri
            waitbar(i/NbTri,wbm);
           A(i,i) = 1/100*sum(Coef_trans(i,:));
           b(i,1) = Volume(1,i)*f(Centre_tri(:,i)');
           for j=1:size(TN,2)
               if isnan(TN(i,j))
                   if Typesegment(1,3*(i-1)+j)==0
                       A(i,i) = A(i,i) + Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1) - 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0);
                   elseif Typesegment(1,3*(i-1)+j)==1
                       A(i,i) = A(i,i) + Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1) - 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0); 
                       b(i,1) = b(i,1) + 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                       b(i,1) = b(i,1) + Coef_trans(i,j)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   elseif Typesegment(1,3*(i-1)+j)==2
                       A(i,i) = A(i,i) - 1/100*Coef_trans(i,j);
                       A(i,i) = A(i,i) + Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1) - 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0); 
                       b(i,1) = b(i,1) + Long_segment(1,3*(i-1)+j)*h(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));%'*([0 -1;1 0]*(Points(Eall(3*(i-1)+j,2),:)-Points(Eall(3*(i-1)+j,1),:))')/norm(Points(Eall(3*(i-1)+j,2),:)-Points(Eall(3*(i-1)+j,1),:),2);
                   end

               else
                   A(i,TN(i,j)) = -1/100*Coef_trans(i,j) + 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0);
                   if Typesegment(1,3*(i-1)+j)==0
                       A(i,i) = A(i,i) + Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1) - 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0);
                   elseif Typesegment(1,3*(i-1)+j)==1
                       A(i,i) = A(i,i) + Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1) - 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0);
                       b(i,1) = b(i,1) + 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                       b(i,1) = b(i,1) + 1/100*Coef_trans(i,j)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   elseif Typesegment(1,3*(i-1)+j)==2
                       A(i,i) = A(i,i) - 1/100*Coef_trans(i,j) ;
                       A(i,i) = A(i,i) + Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1) - 1/Volume(1,i)*Long_segment(1,3*(i-1)+j)*Prodvitnor_seg(3*(i-1)+j,1)*(Prodvitnor_seg(3*(i-1)+j,1)>0);
                       b(i,1) = b(i,1) + Long_segment(1,3*(i-1)+j)*h(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   end
               end
            end
        end
    elseif equation==4
        beta = 5e-2;%input('beta = ');
        D = 1e-4;%input('D = ');
        A = sparse(NbTri,NbTri);
        b = sparse(NbTri,1);
        wbm = waitbar(0,'Building A and b');
        for i=1:NbTri       
            waitbar(i/NbTri,wbm);
           A(i,i) = D*sum(Coef_trans(i,:))+beta*Volume(1,i);
           b(i,1) = Volume(1,i)*f(Centre_tri(:,i)');
           for j=1:size(TN,2)
               if isnan(TN(i,j))
                   if Typesegment(1,3*(i-1)+j)==0

                   elseif Typesegment(1,3*(i-1)+j)==1
                       b(i,1) = b(i,1) + D*Coef_trans(i,j)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   elseif Typesegment(1,3*(i-1)+j)==2
                       A(i,i) = A(i,i) - D*Coef_trans(i,j);
                       b(i,1) = b(i,1) + D*Long_segment(1,3*(i-1)+j)*h(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));%'*([0 -1;1 0]*(Points(Eall(3*(i-1)+j,2),:)-Points(Eall(3*(i-1)+j,1),:))')/norm(Points(Eall(3*(i-1)+j,2),:)-Points(Eall(3*(i-1)+j,1),:),2);
                   end

               else
                   A(i,TN(i,j)) = -D*Coef_trans(i,j);
                   if Typesegment(1,3*(i-1)+j)==0

                   elseif Typesegment(1,3*(i-1)+j)==1
                       b(i,1) = b(i,1) + D*Coef_trans(i,j)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   elseif Typesegment(1,3*(i-1)+j)==2
                       A(i,i) = A(i,i) - D*Coef_trans(i,j) ;
                       b(i,1) = b(i,1) + D*Long_segment(1,3*(i-1)+j)*h(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   end
               end
           end
        end
    elseif equation==5
        A = sparse(NbTri,NbTri);
        b = sparse(NbTri,1);
        D = 5e-2;
        r = 0;
        wbm = waitbar(0,'Building A and b');
        for i=1:NbTri       
           waitbar(i/NbTri,wbm);
           A(i,i) = D*sum(Coef_trans(i,:))-r*Volume(1,i);
           b(i,1) = Volume(1,i)*f(Centre_tri(:,i)');
           for j=1:size(TN,2)
               if isnan(TN(i,j))
                   if Typesegment(1,3*(i-1)+j)==0

                   elseif Typesegment(1,3*(i-1)+j)==1
                       b(i,1) = b(i,1) + D*Coef_trans(i,j)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   elseif Typesegment(1,3*(i-1)+j)==2
                       A(i,i) = A(i,i) - D*Coef_trans(i,j);
                       b(i,1) = b(i,1) + D*Long_segment(1,3*(i-1)+j)*h(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));%'*([0 -1;1 0]*(Points(Eall(3*(i-1)+j,2),:)-Points(Eall(3*(i-1)+j,1),:))')/norm(Points(Eall(3*(i-1)+j,2),:)-Points(Eall(3*(i-1)+j,1),:),2);
                   end

               else
                   A(i,TN(i,j)) = -D*Coef_trans(i,j);
                   if Typesegment(1,3*(i-1)+j)==0

                   elseif Typesegment(1,3*(i-1)+j)==1
                       b(i,1) = b(i,1) + D*Coef_trans(i,j)*g(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   elseif Typesegment(1,3*(i-1)+j)==2
                       A(i,i) = A(i,i) - D*Coef_trans(i,j) ;
                       b(i,1) = b(i,1) + D*Long_segment(1,3*(i-1)+j)*h(0.5*(Points(Eall(3*(i-1)+j,1),:)+Points(Eall(3*(i-1)+j,2),:)));
                   end
               end
           end
        end
        
    end
    close(wbm);
end