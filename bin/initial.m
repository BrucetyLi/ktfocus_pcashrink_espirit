function W = initial(pic,V)
[row,col,frame] = size(pic);
W = reshape(pic,[row*col,frame])*conj(V);
W = reshape(W,size(pic));
end


