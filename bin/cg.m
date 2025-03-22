function x = cg(A,b,x0,steps)
grad = A(x0) -b;
coju = -grad;
step = 0;
eps = 1.0e-6;

while 1 
    if norm(grad(:)) < eps
        break;
    end
    if step >steps
        break;
    end
    step = step +1;
    temp1 = A(coju);
    temp = conj(coju).*temp1;
    alpha = norm(grad(:))^2/sum(abs(temp(:)));
    x1 = x0 + alpha*coju;
    grad1 = grad + alpha*temp1;
    beta = norm(grad1(:))^2/norm(grad(:))^2;
    coju1 = -grad1 + beta*coju;

    grad = grad1;
    coju = coju1;
    x0 = x1;
    %disp(norm(grad(:)));
end
x = x0;
    
end


