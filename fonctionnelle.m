function F = fonctionnelle(x,f,f_prime,g,g_prime,h,h_prime,TN,Coef_trans,Volume,C,Delta_t,U)
    F = zeros(size(TN,1),1);
    for i=1:size(TN,1)
        F(i,1) = -h(x(i))*Delta_t+x(i)-U(i);
        for j=1:size(TN,2)
            if isnan(TN(i,j))
                
            else
                F(i,1) = F(i,1)+(x(i)-x(TN(i,j)))*Coef_trans(i,j)*0.5*(f(x(i))+f(x(TN(i,j))))*Delta_t/Volume(1,i) ...
                    + Coef_trans(i,j)*g(x(i))*(C(TN(i,j))-C(i))*(C(TN(i,j))-C(i)>0)*Delta_t/Volume(1,i)...
                    + Coef_trans(i,j)*g(x(TN(i,j)))*(C(TN(i,j))-C(i))*(C(TN(i,j))-C(i)<0)*Delta_t/Volume(1,i);
            end
        end
    end
end