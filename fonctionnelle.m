function F = fonctionnelle(x,f,f_prime,g,g_prime,h,h_prime,TN,Coef_trans,Volume,C,Delta_t,U)
    F = zeros(size(TN,1),1);
    for i=1:size(TN,1)
    F(i,1) = x(i) - Delta_t*h(x(i)) - U(i);
        for j=1:size(TN,2)
            if isnan(TN(i,j))
                
            else
                F(i,1) = F(i,1) + Delta_t/Volume(1,i)*Coef_trans(i,j)*(f(x(i))-f(x(TN(i,j)))) + Delta_t/Volume(1,i)*Coef_trans(i,j)*g(x(i))*max(C(TN(i,j),1)-C(i,1),0)-Delta_t/Volume(1,i)*Coef_trans(i,j)*g(x(TN(i,j)))*max(C(i,1)-C(TN(i,j),1),0);
            end
        end
    end
end