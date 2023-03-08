function [W,pxf] = demandweight(w,D,in)
% Proximity factor
   % Adjust to make sum_i w_ij = w_j  (inbound weight)
proxfac2 = @(D,w,p) proxfac3(D,w,p).*(w(:)'./sum(proxfac3(D,w,p),1));
fh = @(p,D,w) sum(sum(proxfac2(D,w,p).*D));
pxf = fminsearch(@(p) abs(fh(p,D,w) - in.dinb),1);
W = proxfac2(D,w,pxf);
if any(any(W < 0))
   error('Negative element(s) detected.')
end
if in.dodisp
   fprintf('Proximity factor = %f for %f avg inbound dist.\n',...
      pxf,in.dinb);
end
