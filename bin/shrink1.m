function y=shrink1(x,lambda)
y=sign(x).*max(abs(x)-lambda,0);
end