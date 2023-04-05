%% Loading Workspace
clear
clc
load('Raleigh80.mat');
noLU = 4; % number of loading/unloading docks
noM = 10*5; % number of modules per array
noA = 8*8-4; % number of arrays per level
ur = 0.8; % utilization rate
%% Arrival Rate
Pop = sum(q); %population
lph = Pop*in.rwk/7/8/7/2; % loads per hour
GSlph = lph/noGS; % loads per hour per grocery store
noDC = length(u); % number of DCs
DClph = GSlph*noGS/noDC % DC loads per hour
DClph1 = DClph/ur; % adjusting for 80% utilization
mu = DClph1/noLU; % server utilization
LUtime = 60/mu % loading/unloading time to satisfy demand or min L/U time
%% Number of levels in DC
npph = 2.51; % number of people per household in the US
dem = round(npph*in.rwk); % demand per house hold
op = 7; % order point. 7 packages triggers an order
ot = round((op/dem)*7); % order time. Days between orders.
for i = 1:(op-1)
    at(i) = ot-(ot/op)*i; % average wait time for each package
end
awt = sum(at/(op)) % average wait time
L = DClph1*in.k*16*awt % average packages in DC
NoL = round(L/(noM*noA)) % number of levels required per DC
%% Cost Analysis
% Module Cost
actvh = round(in.dinb/(in.v/60/(DClph/60))); % active vehicles per DC
mod = (noDC*noA*noM*NoL) + in.k*actvh*noDC; % number of modules needed for the entire region
mcost = 50; % module cost, currently set as $50
Totmcost = mcost*mod;
salvperc = 0.25; % salvage percentage
salv = salvperc*Totmcost; % salvage value
Oppcost = Totmcost*1.02^5; % including 2% interest rate over 5 years
Modcost = ((Oppcost - salv)/5)/12 % monthly cost for modules per city
Modcost2 = ((Totmcost - salv)*(0.02/(1-(1+0.02)^-5))+salv*0.02)/12

% DC electricity cost
aweu = 6.6; % average warehouse electricity usage per square foot in kwh, annual
msize = 0.25; % module size
Aw = 0.25^2*noA*noM*NoL*10.76391; % total area of DC in sq ft
Ew = aweu*Aw; % Total annual electricity usage per DC
e = 0.0926; % commercial electricity price in that region in $
ToteDC = (Ew/12)*e; % electricity cost per DC per month
Tote = ToteDC*noDC; % monthly cost for DC electricity usage per city

% Vehicle leasing cost
lcost = 2750; % leasing cost per vehicle per month
Novh = actvh*noDC; % number of vehicles needed for the entire region
Totvh = lcost*Novh; % monthly leasing cost for driverless delivery vehicles

% Vehicle electricity cost
vhe = 280/1000; % vehicle electricity usage kwh/mi
DCpkg = (in.rwk*Pop*30/7)/noDC; % number of packages per DC per month
vhpkg = DCpkg/actvh; % number of packages per vehicle per month
vhld = vhpkg/in.k; % number of loads per vehicle per month
vhmi = vhld*in.dinb; % average distance traveled per vehicle per month (mi)
Totveu = vhmi*vhe; % total monthly electricity usage per vehicle
Totv = Totveu*e*Novh; % monthly total vehicle electricity cost

% Real Estate Cost
rc = 2; % monthly rent per square foot for commercial space
Acw = Aw*3/NoL; % actual warehouse area as we are planning to fit in 3 levels into one
TotrDC = Acw*rc; % monthly total rent per DC
Totr = TotrDC*noDC; % monthly total rent cost for all DCs

Totc = Totr + Totv + Totvh + Tote + Modcost2; % monthly total cost
Totd = (in.rwk*Pop*30/7)/7; % number of loads per month which is more or less number of deliveries
Dcost = Totc/Totd % cost per delivery