function [fP,fPD,fPH,fPS] = initialloaddemand(q,in)
fPD = in.r*q/in.k; % D for delivery ld/hr
fPH = in.r*q/(in.k*2); % H for linehaul ld/hr
fPS = in.r*q/in.k; % S for supply store to a DC k
fP = (fPD + fPH + fPS)'; % f stands for number of loads per hour
vdisp('sum(fP)')