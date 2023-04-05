function [X,A,Aw] = HDLN2030(CSAstr,in)
%% Select census block demand points around metro region
[XP,q,aP,zADP] = getADPforCSA30(CSAstr,in);
%rADP = sqrt(sum(zADP.LandArea)/pi) % Radius of ADP region
noGS = round(sum(q)/in.qmin); % number of grocery stores, 10000 comes from my calculations*
dlb = 0.376*sqrt(sum(aP)/noGS); % d lower bound, we can cite this later
dub = 0.51*sqrt(sum(aP)/noGS); % d upper bound, we can cite this later
in.dinb = sqrt(dlb*dub) % Avg dist (mi) inbound to DC
%%
[f,fPD,fPH,fPs] = initialloaddemand(q,in);
fP = f; P = XP; aP = zADP.LandArea';

WP = demandweight(q/sum(q),dagg(P,P,aP),in); % output of the proximity factor
tic
clear FP uerr
i = 1;
[done,X] = checkfmax2(fP,in);
while ~done
    fprintf('i = %d\n',i)
    FP(i,:) = fP; % current estimate for load per hour
    [X,F,out1] = locateDC(fP,XP,aP,in); % Calculating the DC locations using this function
    l(i) = size(X,1); % storing the number of DC location data
    A = F./fP; % Z matrix shows how much a location is being served by which DC. If there is a split, the value is less than one, if it is served full it will be one and if it is not being served than it will be zero
    W = A*WP*A'; % It collects the inbound loads for each ADP and then allocates the sum to the outbound loads, more like a proximity factor output for DCs
    fL = localloaddemand(q,W,in);
    a = sum(aP(:)'.*A,2);
    if size(X,1) == 1 % try to find a solution for this? incorporate try catch?
        fH = 0;
        Aw = W;
        T = 0;
    else
        [fH,T,Aw] = linehaulloaddemand(X,W,sum(q),a,in);
    end
    f = fL + fH; % recalculated load/hr per DC
    u = f/in.fmax;  % Load capacity factor at DC
    if in.dodisp
      vdisp('fL,fH,f,u',1)
    end
    uerr(i) = sum(.90 - u(u<.90)) + sum(u(u>1.1) - 1.1);
    fprintf('uerr = %f\n',uerr(i))
    if sum(l==1) > 3 && sum(l==2) > 3
    warning('Current demand is too large for one DC and too small for two DCs. Relaxation applied')
        if all(u > 0.9 - i*0.01 & u < 1.1 + i*0.01)
            disp('Feasible soln reached.')
            done = true;
        elseif i > 100
            disp('Max iter reached.')
            done = true;
        else
            i = i+1;
            fP1 = u' * (A.*fP); % i may need to investigate this and the lines after this
            alpha = 0.8; %(sqrt(5) - 1)/2;  % Golden ratio conjugate
            fP = alpha*fP + (1 - alpha)*fP1; % 80 percent of actual total load, 20 percent of the solution of the iteration
            [done,X] = checkfmax2(fP,in);
            sum(fP)
        end
    else
    if all(u > 0.9 & u < 1.1) % that's why we have 1.1 and 0.9, because we don't need to hit exact 100% since we are using golden rule
      disp('Feasible soln reached.')
      done = true;
   elseif i > 100
      disp('Max iter reached.')
      done = true;
   else
       i = i+1;
       fP1 = u' * (A.*fP); % i may need to investigate this and the lines after this
       alpha = 0.8; %(sqrt(5) - 1)/2;  % Golden ratio conjugate
       fP = alpha*fP + (1 - alpha)*fP1; % 80 percent of actual total load, 20 percent of the solution of the iteration
       [done,X] = checkfmax2(fP,in);
       sum(fP)
    end
    end
end
toc
%%
if X ~= 0
if in.dodisp, plotlognet(X,XP,A,Aw,zADP.citystr); end %plotLHlanes(Aw,X); end
end
%%
if X ~= 0
TD = dagg(X,XP,aP(:)')/in.v;
WD = fPD(:)'.*A./sum(fPD);
tD = sum(sum(WD.*TD))*60
tHh = sum(sum((W/(1 - sum(diag(W)))).*T))*60
sumWH = 1 - sum(diag(W))
end