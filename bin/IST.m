function w =  IST(x,T,y,lambda)
w = shrink1(x + T(y-T(x,1),2),lambda/2);

end


