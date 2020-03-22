function gradF = gradient_fonctionnelle(x,f,f_prime,g,g_prime,h,h_prime,TN,Coef_trans,Volume,C,Delta_t)
    gradF = sparse(zeros(size(TN,1),size(TN,1)));
    for i=1:size(TN,1)
        gradF(i,i) = 1 - Delta_t*h_prime(x(i));
        for j=1:size(TN,2)
            if isnan(TN(i,j))
                
            else
                gradF(i,i) = gradF(i,i) + Delta_t/Volume(1,i)*Coef_trans(i,j)*f_prime(x(i)) + Delta_t/Volume(1,i)*Coef_trans(i,j)*g_prime(x(i))*(C(TN(i,j),1)-C(i,1));%max(C(TN(i,j),1)-C(i,1),0) ;
                gradF(i,TN(i,j)) = -Delta_t/Volume(1,i)*Coef_trans(i,j)*f_prime(x(TN(i,j))) - Delta_t/Volume(1,i)*Coef_trans(i,j)*g_prime(x(TN(i,j)))*(C(TN(i,j),1)-C(i,1));%*max(C(i,1)-C(TN(i,j),1),0);
            end
        end
    end
end