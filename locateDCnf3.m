function [X,F,TC,X30new] = locateDCnf3(f,P,a,X1,X30,in)
% This function is for adding new facilities to the existing facilities

% Determine DC locations to serve ADPs.
% [X,Q,A] = getDCloc(wL,wH,P,q,a,in)

K = sum(f)/round(sum(f)/in.fmax); % capacity per DC (this should be very close to fmax) this is done for exact allocation
% D = dists(P,P,'mi'); 
D = dagg(P,P,a); %aggregate distance matrix for census block group locations

% Existing facility portion
D1 = dists(X1,P);
ef = argmin(D1,2);
Fef = zeros(1,length(f));
for i = 1:length(X1)
    idx11 = argsort(D(ef(i),:)); % index ranking according to distances for the row of starting point
    idx11(is0(f(idx11))) = []; % deleting zero load values in f matrix
    csum1 = cumsum(f(idx11)); % cumulative sum of the loads of datapoints that are part of idx vector
    idx22 = find(csum1 <= K,1,'last'); % find the point where we cross K~96 indicating that we reached the DC capacity.
    Fef(i,idx11(1:idx22)) = f(idx11(1:idx22)); % store the served load into F matrix
    f(idx11(1:idx22)) = 0; % delete the already served census block locations
    if ~is0(csum1(idx22) - K)
      frem1 = K - csum1(idx22);  % f of last EF alloc to NF
      idx22 = idx22 + 1;
      Fef(i,idx11(idx22)) = frem1; % allocating the uncapped capacity to the next closest point
      f(idx11(idx22)) = f(idx11(idx22)) - frem1; % taking out the served load from the data point
    end
    Xef(i,:) = X1(i,:); % finding the coordinates of the DC
end

% F function shows which locations are served by the corresponding DC and the load that has been served. The
% row number also incdicates the DC number.

% X function shows the coordinates of DCs

% TC0 shows the lowest TC

[X,F,TC] = getDCloc1(f,P,a,D,K);
D30 = dists(X30,X,'mi');
[M,I] = min(D30);
i = 1;

while length(I)~=length(unique(I))
    [X,F,TC] = getDCloc1(f,P,a,D,K);
    D30 = dists(X30,X,'mi');
    [M,I] = min(D30);
    i = i + 1;
    if i == 1000 % if we try 1000 times and still can't allocate all DCs separately we manually separate them
        for k = 1:size(D30,2)
            [M,I1] = min(D30(:,k));
            I(1,k) = I1;
            D30(I1,:) = Inf;
        end
    end
end

XX = X30(I,:);
X30(I,:) = [];
X30new = X30;

X = [XX; Xef];
F = [F;Fef];

%%%%%%%%
function [X,F,TC] = getDCloc1(f,P,a,D,K)
% disp([length(f)])


C = D.*f; % cost in terms of load times aggregate distance
idxCH = convhull(P); idxCH(end) = []; % this would give us the index values of convex hull
% idxCHi = idxCH(1);
idxCHi = idxCH(wtrandperm(f(idxCH),1)); % wtrandperm selects initial point on CH randomly in proportion to its load demand, gives different results in each run

F = zeros(1,length(f));
X = [0 0];
i = 1;

TC = 0;
done = false;
while ~done
   idx = argsort(D(idxCHi,:)); % index ranking according to distances for the row of starting point
   idx(is0(f(idx))) = []; % deleting zero load values in f matrix
   csum = cumsum(f(idx)); % cumulative sum of the loads of datapoints that are part of idx vector
   idx2 = find(csum <= K,1,'last'); % find the point where we cross K~96 indicating that we reached the DC capacity.
   F(i,idx(1:idx2)) = f(idx(1:idx2)); % store the served load into F matrix
   f(idx(1:idx2)) = 0; % delete the already served census block locations
   if ~is0(csum(idx2) - K)
      frem = K - csum(idx2);  % f of last EF alloc to NF
      idx2 = idx2 + 1;
      F(i,idx(idx2)) = frem; % allocating the uncapped capacity to the next closest point
      f(idx(idx2)) = f(idx(idx2)) - frem; % taking out the served load from the data point
   end
   [TCi,idx3] = min(sum(C(idx(1:idx2),idx(1:idx2)),2)); % from the locations that are served, we are picking the lowest load*distance value for our DC
   X(i,:) = P(idx(idx3),:); % finding the coordinates of the DC
   TC = TC + TCi; % total cost of DC added
   
   idxCH(is0(f(idxCH))) = []; % if some of the convex hull data points are fully served, we get rid of them
   if ~isempty(idxCH) % checking if we still have unserved convex hull perimeter points
      idxCHi = idxCH(argmin(D(idxCHi,idxCH))); % taking the shortest distance from our ith iteration as the starting point of i+1th iteration
      i = i + 1;
   else
      if any(~is0(f)) % checking whether there are still unserved locations or not
         is = ~is0(f); % index vector for unserved locations
         % THIS PART MIGHT BE SLOWING OUR CODE DOWN
         [X1,F1,TC1] = getDCloc1(f(is),P(is,:),a(is),D(is,is),K); % rerun the function for the other unserved locations again.
         X = [X; X1]; % save the calculated data into the existing vectors
         F(end+1:end+size(F1,1),is) = F1;
         TC = TC + TC1;
      end
      done = true; % if all the data points are served, exit the while loop
   end
end