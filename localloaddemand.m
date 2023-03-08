function fL = localloaddemand(q,W,in)
fD = in.r*sum(q)*sum(W,1)'/in.k;
fS = in.r*sum(q)*sum(W,2)/in.k;
fL = fD + fS;