function [P,q,a,x] = getADPforCSA(CSAstr,in)
% Select census block demand points around metro region.
% x = getADPforCSA(CSAstr,in)

[SCfips,Name,str] = uscounty('SCfips','Name','CSAtitle',...
   mor({CSAstr},uscounty('CSAtitle')));
str = [str{1} ' CSA'];
if in.dodisp, mdisp(SCfips',Name,{'SCfips'},str), end
x = uscenblkgrp(mor(SCfips,uscenblkgrp('SCfips')));
x = subsetstruct(x,x.Pop ~= 0);
x.Pop(x.Pop>in.qmaxADP) = in.qmaxADP;  % pop > qmax => fP > fmax for large data points
xycog = @(XY,w) [w(:)'*XY(:,1) w(:)'*XY(:,2)]/sum(w); % center of gravity
xyCOG = xycog(x.XY,x.Pop);
s = x.Pop./x.LandArea; % density
if in.dodisp, vdisp('min(s),max(s),mean(s),median(s)'), end

iss = s > 500; % indexing for high density data points / using 500 like we did in tmax calculation
%
dt = delaunayTriangulation(x.XY); % or voronoi diagram
IJ = abs(tri2list(dt.ConnectivityList));
IJ(~iss(IJ(:,1)) | ~iss(IJ(:,2)),:) = []; % we get rid of the low density points
[t,nk] = minspan([IJ ones(size(IJ,1),1)]);
dmin = Inf;
for i = 1:nk % maybe put an if statement before this to enhance the code?
   idxi = IJ(t == i,:); idxi = unique(idxi(:));
   di = dists(xyCOG,xycog(x.XY(idxi,:),x.Pop(idxi)),'mi'); % we are selecting the tree that is closest to the COG
   if di < dmin, dmin = di; idxmin = idxi; end
end
%
k = convhull(x.XY(idxmin,:));
idx = find(inpolygon(x.XY(:,1),x.XY(:,2),...
   x.XY(idxmin(k),1),x.XY(idxmin(k),2))); % find the index values of points inside the convex hull
is = inpolygon(sub(uscity('XY'),'#(:,1)'),...
   sub(uscity('XY'),'#(:,2)'),x.XY(idx(k),1),x.XY(idx(k),2));
if any(is)  % empty for Charlotte
s = uscity(is);
% citystr = s.Name(argmax(s.Pop)); tmp = s.ST(argmax(s.Pop));
% citystr = [citystr{:} ', ' tmp{:}]
citystr = [sub(s.Name(argmax(s.Pop)),'#{:}') ', ' ...
   sub(s.ST(argmax(s.Pop)),'#{:}')];
else
   citystr = CSAstr;
end
if in.dodisp
   makemap(x.XY),pplotroads(x.XY)
   k = convhull(x.XY(idx,:));
   h1 = pplot(x.XY(iss,:),'r.','DisplayName','High Density ADP');
   h2 = pplot(x.XY(~iss,:),'k.','DisplayName','Low Density ADP');
   h3 = pplot(x.XY(idx(k),:),'b-','LineWidth',1,'DisplayName',citystr);
   title(str)
   legend([h1 h2 h3])
   vdisp('sum(x.Pop),sum(x.Pop(idx)),sum(x.Pop(idx))/sum(x.Pop)')
end
x = subsetstruct(x,idx);
s = x.Pop./x.LandArea;
if in.dodisp, vdisp('min(s),max(s),mean(s),median(s)'), end
x.CSAstr = str;
x.citystr = citystr;
[P,q,a] = deal(x.XY,x.Pop,x.LandArea);
