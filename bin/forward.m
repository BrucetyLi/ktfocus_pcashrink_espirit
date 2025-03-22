
% function Ax = forward(x,U,opt)
% x = reshape(x,size(U));
% if opt ==1
% Ax = U.*fftc(ifftc(x,3),1);
% else
% Ax = ifftc(fftc(U.*x,3),1);
% end
% Ax = Ax(:);
% end


function Ax = forward(x,U,V,S,opt)
[row,col,frame] = size(U);
channel = size(S,3);
shape = [row,col,frame,channel];
if opt == 1
 temp =repmat(reshape(reshape(x,[row*col,frame])*conj(V'),size(U)),[1,1,1,channel]);
 temp1 = temp.*repmat(reshape(S,[row,col,1,channel]),[1,1,frame,1]);
Ax = fftc(temp1,1).* repmat(U,[1,1,1,channel]);
else
 x = reshape(x,shape);
 temp2 = repmat(reshape(conj(S),[row,col,1,channel]),[1,1,frame,1]).*ifftc(repmat(U,[1,1,1,channel]).*x,1);
Ax = reshape(sum(temp2,4),[row*col,frame])*conj(V);
end
Ax = Ax(:);
end


