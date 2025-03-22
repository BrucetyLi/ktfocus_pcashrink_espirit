function b = right(W,v)
b = W.*ifftc(fftc(v,3),1);
b = b(:);
end


% function b = right(C,U,v,p)
% [row,col,frame,channel] = size(C);
% b = zeros(row,col,frame,channel);
% for j = 1:col
% for i = 1:channel
%     C1 = squeeze(C(:,j,:,i));
%     U1 = squeeze(U(:,j,:,i));
%     b(:,j,:,i) = U1.*fftc(ifftc(C1.*squeeze(p(:,j,:)),2),1);
% 
% end
% end
% b  =-b(:)+v(:);
% end




