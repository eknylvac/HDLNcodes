function [fH,T,Aw] = linehaulloaddemand(X,W,Q,a,in)

fHh = @(Aw) in.r*Q*(sum(Aw,1)' + sum(Aw,2))/(in.k*2); % summing up the W matrix without the diagonal
tHh = @(W,T) sum(sum((W/(1 - sum(diag(W)))).*T)); % LH transit (hr) incl LU

D = dagg(X,X,a); % aggregate distance between DCs, here the area is allocated depending upon how much is served by which DC

At = (triu(D,1) + tril(D,-1))/in.v; % Linehaul arc travel time (hr), diagonal is zero
Aw = triu(W,1) + tril(W,-1); % The W is to remove set the local-only demand to 0 so that it is not included in the linehaul calculations
fH = fHh(Aw); % calculating the linehaul load/hr
A = getarctime(Aw,At,fH,in); % arc time
tH = tHh(W,A);  % T = A for complete graph
done = false;
while ~done
   [T1,Av1] = getLHtime(A,W,in);
   tH1 = sum(sum(W.*T1));
   if tH1 < tH
      tH = tH1; T = T1; Aw = Av1;
      fH = full(fHh(Aw));
      A = getarctime(Aw,At,fH,in);
   else
      done = true;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,out] = getarctime(Av,At,fH,in)
% A = getarctime(Av,At,Q,in)

A = sparse(length(Av),length(Av),0);
A(Av~=0) = At(Av~=0); % arc travel time matrix where Aw is zero A is also zero
Ad = Av;
Ad(Av~=0) = 1./(Av(Av~=0) * sum(fH)); % hr/load, calculating the headway
Ad = Ad * sqrt(0.5*5.5);  % Wait for truck delay (headway delay time)
A = A + Ad;
out.Ad = Ad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,Av,out] = getLHtime(A,W,in)

[T,P] = dijk(A);  % T = travel + L/U time + 
%T = triu(T) + triu(T,1)';  % To make symmetric
P = pred2path(P);
Av = sparse(length(W),length(W),0); %Aw = all paths total weight adj matrix
for i = 1:length(Av)
   for j = 1:length(Av)
%       if i > j, P{i,j} = fliplr(P{j,i}); end  % To make symmetric
      p = P{i,j};  % When i == j, 0 added to Av
      Av = Av + sparse(p(1:end-1),p(2:end),W(i,j),length(Av),length(Av));
   end
end
out.P = P;
