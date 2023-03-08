function [X,F,TC,TC0] = locateDC(f,P,a,in)
% Determine DC locations to serve ADPs.
% [X,Q,A] = getDCloc(wL,wH,P,q,a,in)

K = sum(f)/round(sum(f)/in.fmax); % capacity per DC (this should be very close to fmax) this is done for exact allocation
% D = dists(P,P,'mi'); 
D = dagg(P,P,a); %aggregate distance matrix for census block group locations

TC0 = Inf;
for i = 1:12 % just a random number
   [Xi,Fi,TCi] = getDCloc1(f,P,a,D,K); % getDCloc1 finds the total cost and location of DCs for the entire system
   if TCi < TC0 % using the lowest total cost run for our output
      X = Xi; F = Fi; TC0 = TCi; % storing the lowest TC found
   end
end
% F function shows which locations are served by the corresponding DC and the load that has been served. The
% row number also incdicates the DC number.

% X function shows the coordinates of DCs

% TC0 shows the lowest TC

TC = 0;
for i = 1:size(X,1) % the number of DCs located by getDCloc1
   is = ~is0(F(i,:)); % is index vector shows which location is served by the corresponding DC.
   TCh = @(xy) dagg(xy,P(is,:),a(is)) * F(i,is)'; % total cost in terms of load*distance?
%    TCh = @(xy) dists(xy,P(is,:),'mi') * F(i,is)';
   [X(i,:),TCi] = fminsearch(TCh,X(i,:)); % finds the lowest total cost for corresponding DC and the minimal cost location for DC
   TC = TC + TCi; % add up the total cost to find the ultimate total cost for our design
end

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