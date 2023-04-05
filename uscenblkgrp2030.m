function xx = uscenblkgrp2030(CSAstr)

[SCfips,Name,str] = uscounty('SCfips','Name','CSAtitle',...
   mor({CSAstr},uscounty('CSAtitle')));
str = [str{1} ' CSA'];
x00 = uscenblkgrp2000(mor(SCfips,uscenblkgrp2000('SCfips')));
x00 = subsetstruct(x00,x00.Pop ~= 0);
x10 = uscenblkgrp(mor(SCfips,uscenblkgrp('SCfips')));
x10 = subsetstruct(x10,x10.Pop ~= 0);
x20 = uscenblkgrp2020(mor(SCfips,uscenblkgrp2020('SCfips')));
x20 = subsetstruct(x20,x20.Pop ~= 0);
%% Finding small group borders
xymax00 = [max(x00.XY(:,1))+0.0001 max(x00.XY(:,2))+0.0001]; % finding area extremes
xymin00 = [min(x00.XY(:,1))-0.0001 min(x00.XY(:,2))-0.0001];

xymax10 = [max(x10.XY(:,1))+0.0001 max(x10.XY(:,2))+0.0001];
xymin10 = [min(x10.XY(:,1))-0.0001 min(x10.XY(:,2))-0.0001];

xymax20 = [max(x20.XY(:,1))+0.0001 max(x20.XY(:,2))+0.0001];
xymin20 = [min(x20.XY(:,1))-0.0001 min(x20.XY(:,2))-0.0001];

incr00 = [(xymax00(1,1)-xymin00(1,1))/10 (xymax00(1,2)-xymin00(1,2))/10];

incr10 = [(xymax10(1,1)-xymin10(1,1))/10 (xymax10(1,2)-xymin10(1,2))/10];

incr20 = [(xymax20(1,1)-xymin20(1,1))/10 (xymax20(1,2)-xymin20(1,2))/10];

xycog = @(XY,w) [w(:)'*XY(:,1) w(:)'*XY(:,2)]/sum(w); % center of gravity

%% Calculation of 2030 population
TotPop00 = zeros(10,10);
TotPop10 = zeros(10,10);
TotPop20 = zeros(10,10);
TotPop30 = zeros(10,10);
for i=1:10
    for j=1:10
        xy00ext = [xymin00(1,1)+i*incr00(1,1) xymin00(1,2)+j*incr00(1,2)
            xymin00(1,1)+(i-1)*incr00(1,1) xymin00(1,2)+(j-1)*incr00(1,2)]; %extreme points of the small group
        idx = find(xy00ext(2,1) < x00.XY(:,1) & x00.XY(:,1) < xy00ext(1,1) & ...
            xy00ext(2,2) < x00.XY(:,2) & x00.XY(:,2) < xy00ext(1,2));
        
        xy00{i,j} = x00.XY(idx,:);
        Pop00{i,j} = x00.Pop(idx,:);
        a00{i,j} = x00.LandArea(idx,:);
        d00{i,j} = x00.Pop(idx,:)./x00.LandArea(idx,:); % density
        if xy00{i,j} ~= 0
            COG00{i,j} = xycog(xy00{i,j},Pop00{i,j}); % center of gravity for all small groups
        else
            COG00{i,j} = [mean(xy00ext(:,1)) mean(xy00ext(:,2))];
        end

        TotPop00(i,j) = sum(Pop00{i,j});

        xy10ext = [xymin10(1,1)+i*incr10(1,1) xymin10(1,2)+j*incr10(1,2)
            xymin10(1,1)+(i-1)*incr10(1,1) xymin10(1,2)+(j-1)*incr10(1,2)]; %extreme points of the small group
        idx = find(xy10ext(2,1) < x10.XY(:,1) & x10.XY(:,1) < xy10ext(1,1) & ...
            xy10ext(2,2) < x10.XY(:,2) & x10.XY(:,2) < xy10ext(1,2));
        
        xy10{i,j} = x10.XY(idx,:);
        Pop10{i,j} = x10.Pop(idx,:);
        a10{i,j} = x10.LandArea(idx,:);
        d10{i,j} = x10.Pop(idx,:)./x10.LandArea(idx,:); % density
        if xy10{i,j} ~= 0
            COG10{i,j} = xycog(xy10{i,j},Pop10{i,j}); % center of gravity for all small groups
        else
            COG10{i,j} = [mean(xy10ext(:,1)) mean(xy10ext(:,2))];
        end

        TotPop10(i,j) = sum(Pop10{i,j});

        xy20ext = [xymin20(1,1)+i*incr20(1,1) xymin20(1,2)+j*incr20(1,2)
            xymin20(1,1)+(i-1)*incr20(1,1) xymin20(1,2)+(j-1)*incr20(1,2)]; %extreme points of the small group
        idx = find(xy20ext(2,1) < x20.XY(:,1) & x20.XY(:,1) < xy20ext(1,1) & ...
            xy20ext(2,2) < x20.XY(:,2) & x20.XY(:,2) < xy20ext(1,2));
        
        xy20{i,j} = x20.XY(idx,:);
        Pop20{i,j} = x20.Pop(idx,:);
        a20{i,j} = x20.LandArea(idx,:);
        d20{i,j} = x20.Pop(idx,:)./x20.LandArea(idx,:); % density
        if xy20{i,j} ~= 0
            COG20{i,j} = xycog(xy20{i,j},Pop20{i,j}); % center of gravity for all small groups
        else
            COG20{i,j} = [mean(xy20ext(:,1)) mean(xy20ext(:,2))];
        end

        TotPop20(i,j) = sum(Pop20{i,j});

        x = [2000 2010 2020];
        y = [TotPop00(i,j) TotPop10(i,j) TotPop20(i,j)];
        yCOG1 = [COG00{i,j}(1) COG10{i,j}(1) COG20{i,j}(1)];
        yCOG2 = [COG00{i,j}(2) COG10{i,j}(2) COG20{i,j}(2)];
        TotPop30(i,j) = interp1(x,y,2030,'linear','extrap');
        COG30{i,j}(1) = interp1(x,yCOG1,2030,'linear','extrap');
        COG30{i,j}(2) = interp1(x,yCOG2,2030,'linear','extrap');
        if TotPop30(i,j) < 0
            TotPop30(i,j) = 0;
        end
    end
end
%% Allocating population change into existing CBGs
PC = TotPop30 - TotPop20; % population change
for i = 1:10 % for loop to allocate population to the existing CBGs for 2030
    for j = 1:10
        if TotPop30(i,j) > 0 % if population value is higher than zero
            mean20 = mean(Pop20{i,j});
            mean30 = TotPop30(i,j)./length(Pop20{i,j});
            diff = round(mean30 - mean20); % we take the mean difference between 30 and 20
            if PC(i,j) <= 0 % if population decreases (diff is a negative value)
                for k = 1:length(Pop20{i,j})
                    chng = diff + round(std(Pop20{i,j}) * randi([-1000 1000])/1000); % expected population change per CBG
                    Pop30{i,j}(k) = Pop20{i,j}(k) + chng; % subtract the change
                    PC(i,j) = PC(i,j) - chng; % add the changes to the population change matrix
                end
                if PC(i,j) < 0
                    add = round(PC(i,j)./length(Pop20{i,j}));
                    for k = 1:length(Pop20{i,j})
                        Pop30{i,j}(k) = Pop30{i,j}(k) + add; % trying to make sure we don't have negative population change left
                        PC(i,j) = PC(i,j) - add;
                    end
                    if PC(i,j) ~= 0 && abs(PC(i,j)) < 100
                        PC(i,j) = 0;
                    end
                end
            else % if population increases (if diff is a positive value)
                for k = 1:length(Pop20{i,j})
                    if Pop20{i,j}(k) > 3000 % if the 2020 population is larger than 3000 we do nothing
                        Pop30{i,j}(k) = Pop20{i,j}(k);
                    else
                        if Pop20{i,j}(k) + diff > 3000 % if the 2030 population will be above 3000 we cap the mean at 3000
                            Pop30{i,j}(k) = 3000 + round(std(Pop20{i,j}) * randi([-1000 1000])/1000);
                            PC(i,j) = PC(i,j) - (Pop30{i,j}(k) - Pop20{i,j}(k));
                        else
                            chng = diff + round(std(Pop20{i,j}) * randi([-1000 1000])/1000); % expected population change per CBG
                            Pop30{i,j}(k) = Pop20{i,j}(k) + chng; % subtract the change
                            PC(i,j) = PC(i,j) - chng; % subtract the changes to the population change matrix
                        end
                    end
                end
                if PC(i,j) < 0 % making sure that there are no negative changes so we reallocate the negative allocations
                    add = round(PC(i,j)./length(Pop20{i,j}));
                    for k = 1:length(Pop20{i,j})
                        Pop30{i,j}(k) = Pop30{i,j}(k) + add; % trying to make sure we allocate all the population change
                        PC(i,j) = PC(i,j) - add;
                    end
                    if PC(i,j) ~= 0 && abs(PC(i,j)) < 100
                        PC(i,j) = 0;
                    end
                end
            end
        else
            Pop30{i,j} = [];
            PC(i,j) = 0;
        end
    end
end
xy30 = xy20; % taking 2020 locations for the starting point
%% Creating new CBGs
for i = 1:10
    for j = 1:10
        a = 1;
        while PC(i,j) ~= 0
            if PC(i,j) < 600
                newCBG = PC(i,j); % if the population is below 600 we just get one new BCG
            elseif PC(i,j) >= 600 && PC(i,j) < 3000
                newCBG = randi([600 PC(i,j)]); % of the excess population is less than 3000 we randomly assign population values
            else
                newCBG = randi([600 3000]); % new BCGs cannot be larger than 3000
            end
            k = length(Pop30{i,j});
            Pop30{i,j}(k+1) = newCBG; % adding population values to the new CBGs
            a = a + 1;
            PC(i,j) = PC(i,j) - newCBG;
        end
    end
end
for i = 1:10 % for loops to find the location of new CBGs
    for j = 1:10
        noPop30 = length(Pop30{i,j});
        noPop20 = length(Pop20{i,j});
        n = length(Pop30{i,j});
        k = length(Pop20{i,j});
        xy20ext = [xymin20(1,1)+i*incr20(1,1) xymin20(1,2)+j*incr20(1,2)
            xymin20(1,1)+(i-1)*incr20(1,1) xymin20(1,2)+(j-1)*incr20(1,2)]; %extreme points of the small group
        if noPop30 - noPop20 > 0
            COG30input = COG30{i,j};
            Pop30input = Pop30{i,j};
            xy20input = xy20{i,j};
            dCOG1 = @(xy01) abs(COG30input(1) - (Pop30input(1:k) * xy20input(1:k,1) + Pop30input(k+1:n) * xy01)/sum(Pop30input));
            dCOG2 = @(xy02) abs(COG30input(2) - (Pop30input(1:k) * xy20input(1:k,2) + Pop30input(k+1:n) * xy02)/sum(Pop30input)); % function to optimize new locations
            xy01 = repmat(mean(xy20{i,j}(:,1)),n-k,1);
            xy02 = repmat(mean(xy20{i,j}(:,2)),n-k,1); % starting point of new CBGs
            lb01 = repmat(xy20ext(2,1),n-k,1);
            ub01 = repmat(xy20ext(1,1),n-k,1);
            lb02 = repmat(xy20ext(2,2),n-k,1);
            ub02 = repmat(xy20ext(1,2),n-k,1);
            xy30new1 = fmincon(dCOG1,xy01,[],[],[],[],lb01,ub01); % trying fmincon to see if this optimization works better than fminsearch
            xy30new2 = fmincon(dCOG2,xy02,[],[],[],[],lb02,ub02); % optimization to find new coordinate locations for new CBGs
            % if fmincon doesn't work i should try delaunay triangulation
            xy30new = [xy30new1 xy30new2];
            xy30{i,j} = [xy30{i,j}; xy30new]; % adding new locations to the existing coordinates
        end
    end
end
a30 = a20; % getting the initial values from 2020 census for existing CBGs
aa = 1;
for i = 1:10
    for j = 1:10
        if TotPop20(i,j) > 0
            for k = 1:length(Pop20{i,j})
                a20long(aa) = a20{i,j}(k);
                Pop20long(aa) = Pop20{i,j}(k); % creating vectors of population and area to use in interpolation
                aa = aa + 1;
            end
        end
    end
end
[Pop20longu,idxa] = unique(Pop20long); % getting rid of the repeating values for the interpolation
a20longu = a20long(idxa);
for i = 1:10
    for j = 1:10
        n = length(Pop30{i,j});
        k = length(Pop20{i,j});
        if TotPop30(i,j) > 0
            for l = 1:(n-k)
                a30{i,j}(end+1) = interp1(Pop20longu,a20longu,Pop30{i,j}(k+l),'linear','extrap'); % interpolating new area values from the 2020 census
            end
        else
            a30{i,j} = [];
        end
    end
end
%% Saving the output
xx.Pop = [];
for i = 1:10
    for j = 1:10
        for k = 1:length(Pop30{i,j})
            if isempty(xx.Pop) == 1 && isempty(Pop30{i,j}) ~= 1
                xx.Pop(1) = Pop30{i,j}(k);
                xx.XY(1,:) = xy30{i,j}(k,:);
                xx.LandArea(1) = a30{i,j}(k);
            else
                xx.Pop(end+1) = Pop30{i,j}(k);
                xx.XY(end+1,:) = xy30{i,j}(k,:);
                xx.LandArea(end+1) = a30{i,j}(k);
            end
        end
    end
end
xx.Pop = xx.Pop';
xx.Pop(xx.Pop < 0) = 0; % there is a negligible number of negative population points. I am making them zero.
xx.LandArea = xx.LandArea';