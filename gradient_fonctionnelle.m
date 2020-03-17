function gradF = gradient_fonctionnelle(x,f,f_prime,g,g_prime,h,h_prime,TN,Coef_trans,Volume,C,Delta_t)
    gradF = sparse(zeros(size(TN,1),size(TN,1)));
    for i=1:size(TN,1)
        gradF(i,i) = -h_prime(x(i))*Delta_t+1;
        for j=1:size(TN,2)
            if isnan(TN(i,j))
                
            else
                gradF(i,i) = gradF(i,i) + Coef_trans(i,j)*0.5*f_prime(x(i))*x(i)*Delta_t/Volume(1,i)...
                    +  Coef_trans(i,j)*g_prime(x(i))*(C(TN(i,j))-C(i))*(C(TN(i,j))-C(i)>0)*Delta_t/Volume(1,i)...
                    + Coef_trans(i,j)*0.5*(f(x(i))+f(x(TN(i,j))))*Delta_t/Volume(1,i)...
                    - Coef_trans(i,j)*0.5*f_prime(x(i))*x(TN(i,j))*Delta_t/Volume(1,i);
                gradF(i,TN(i,j)) = Coef_trans(i,j)*0.5*f_prime(x(i))*x(i)*Delta_t/Volume(1,i)...
                    - Coef_trans(i,j)*0.5*(f(x(i))+f(x(TN(i,j))))*Delta_t/Volume(1,i)...
                    - Coef_trans(i,j)*0.5*f_prime(x(TN(i,j)))*x(TN(i,j))*Delta_t/Volume(1,i) ...
                    - Coef_trans(i,j)*g_prime(x(TN(i,j)))*(C(TN(i,j))-C(i))*(C(TN(i,j))-C(i)<0)*Delta_t/Volume(1,i);
            end
        end
    end
end