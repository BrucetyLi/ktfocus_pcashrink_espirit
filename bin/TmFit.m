function [T2, B1] = TmFit(S,I,opt)

    [X,~,~,~,~] = lsqnonlin(@diff_sig,opt.lsq.Icomp.X0([1,3]),...
        opt.lsq.Icomp.XL([1,3]),opt.lsq.Icomp.XU([1,3]),opt.lsq.fopt);
    T2 = X(1);  B1 = X(2);
      
function delta = diff_sig(X)
    
     Sfit = I*FSEsig(X(1),X(2),opt);
     delta = abs(S(:) - Sfit(:));
    
end
end

