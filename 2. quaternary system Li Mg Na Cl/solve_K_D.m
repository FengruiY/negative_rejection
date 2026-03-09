function [Di_p,Ki_c] = solve_K_D(lambda,Difu,n)

for i=1:n
    Ki_c(i) = (1+3.867*lambda(i)-1.907*lambda(i)^2-0.834*lambda(i)^3)/(1+1.867*lambda(i)-0.741*lambda(i)^2);

    if lambda(i)<=0.95
        Ki_d(i) = (1.0 + (9/8)*lambda(i)*log(lambda(i)) - 1.56034*lambda(i) + 0.528155*lambda(i)^2 + 1.91521*lambda(i)^3 -2.81903*lambda(i)^4 + 0.270788*lambda(i)^5 + 1.10115*lambda(i)^6 - 0.435933*lambda(i)^7)/(1-lambda(i))^2;
    else
        Ki_d(i) = 0.984*((1-lambda(i))/lambda(i))^(5/2);
    end

    Di_p(i) = Ki_d(i)*Difu(i);
end

end