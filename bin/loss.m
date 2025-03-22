function err = loss(y,T,x,lambda)
err = norm(y- T(x,1)).^2  + lambda*sum(abs(x));
end