% Benchmark with TOUGH3
%% Author : Harun Rashid (hrashid.lsu@gmail.com)
%% Add necessary MRST modules
clc;clear;close all;

mrstModule add ad-props ad-core ad-blackoil compositional geothermal mrst-gui

%% Set up the grid
physdim  = [240 200,0.04];              % domain size in x,y directions 
celldim  = [12,10,1];               % number of cells in x, y dims 
G = cartGrid(celldim, physdim);    % Cartesian Grid 
G = computeGeometry(G);            % grid geometry and connectivity 

%% Find injection faces and monitoring cell indices
% monitoring point near injection
Cellinj = find(G.cells.centroids(:,1)==10 & G.cells.centroids(:,2)==50);   
% Monitoring point in the center
Cellprod = find(G.cells.centroids(:,1)==230 & G.cells.centroids(:,2)==150); 

% monitor=find(G.cells.centroids(:,1)==230 & G.cells.centroids(:,2)==10); 
monitor=60;
G.cells.well=false(G.cells.num,1);
G.cells.well([Cellinj;Cellprod])=true;
% 
% G.cells.wellmntr=false(G.cells.num,1);
% G.cells.wellmntr(Mcen)=true;
% wellID=[Cellinj;Cellprod];
% wellID([Cellinj;Cellprod])=1;
% figure,
% plotGrid(G,'faceColor','y');
% hold on 
% plotGrid(G,wellID,'faceAlpha',1,'faceColor','b');

%% Set fluid structure properties
rhoWS = 1000;
pres=98*barsa;
% Define fluid structure
fluid = initSimpleADIFluid('mu'    , 1.0e-3, ...
    'rho'   , rhoWS , ...
    'phases', 'W',...
    'pRef',    pres);
% pRef = pres;
% c_w = 5e-5/barsa;
% c_o = 1e-5/barsa;
% c_g = 1e-3/barsa;
% 
% fluid.bW = @(p) exp((p - pRef)*c_w);
                       
fluid = addThermalFluidProps(fluid           , ... % Original fluid
                             'Cp'     , 4.2e3, ... % Specific heat capacity
                             'lambdaF', 0.6  , ... % Thermal conductivity
                             'useEOS' , true );    % Use equation of state

%% Make rock structure
perm = 200e-14;
poro = 0.5;
% define rock structure
rock = makeRock(G, perm, poro);
% add thermal props
rock = addThermalRockProps(rock           , ... % Original rock
                           'CpR'    , 1000, ... % Specific heat capacity
                           'lambdaR', 0   , ... % Thermal conductivity
                           'rhoR'   , 2650, ... % Rock density
                           'tau'    , 1   );    % Tortuosity



% Make model
model = GeothermalModel(G, rock, fluid);
model.extraStateOutput = true;
model.outputFluxes     = true;

% Initial state
% HR Edit: Find the temperature using geothermal gradient
Tres = 573.15; %473  temp at top 50k/km geothermal gradient

Tinj = 295.15; %60F ambient T
state0   = initResSol(G, pres, 1);   %(G, 300*barsa, 1);
% state0.T = Tres;
state0.T = ones(G.cells.num,1).*Tres;       % temperature
%% Set up wells
time =1.58E+06;
nstep=35;
ToughTime=[1.00E+02
3.00E+02
7.00E+02
1.50E+03
3.10E+03
6.30E+03
1.27E+04
2.55E+04
5.11E+04
7.67E+04
1.28E+05
1.79E+05
2.82E+05
3.84E+05
4.86E+05
6.91E+05
8.96E+05
1.10E+06
1.31E+06
1.51E+06
1.58E+06
];

dt = diff([0;ToughTime]);
enthalphy = model.fluid.hW(pres,Tinj);
rw=0.5; % Well radius
Inj=Cellinj;
prod1=Cellprod;

% src = addSource([], Inj, 3.5e-3);
% src = addSource(src, prod1, -3.5e-3);
% src = addThermalSource(src, Inj, Tinj);

W = addWell([], G, rock, Inj, ...
    'type', 'rate', 'val',1.e-3, 'sign',1,'compi', 1, 'name', 'inj','Radius', rw,'Dir','z');
% W = addWell([], G, rock, Inj, ...
%     'type', 'bhp', 'val',100*barsa, 'sign',1,'compi', 1, 'name', 'inj','Radius', rw,'Dir','z');
W = addWell(W, G, rock, prod1, ...
    'type', 'bhp', 'val',96.5*barsa,'WI',4.E-12,'sign',-1, 'compi', 1, 'name', 'prod_1','Radius', rw,'Dir','z');


W = addThermalWellProps(W, 'T', Tinj);
% Make schedule
schedule = simpleSchedule(dt, 'W', W);
% schedule = simpleSchedule(rampupTimesteps(time, time/nstep), 'W', W);
%% Run simulation
gravity reset on;
[wellSols, states] = simulateScheduleAD(state0, model, schedule);
% [~, states] = simulateScheduleAD(state0, model, schedule);


figure, 
plotToolbar(G, states)
view(40,30);
axis tight equal;

figure,
plotWellSols(wellSols,cumsum(schedule.step.val))


figure, 
plotToolbar(G, states)
view(40,30);
axis tight equal;
timeSeries=cumsum(schedule.step.val);
for i=1:length(timeSeries)
    monitorpressure(i)=states{i, 1}.pressure(prod1);
end
figure,
plot(timeSeries,monitorpressure,'bo','LineWidth',2)
legend('MRST')
ylabel('Pressure, Pa')
xlabel('Time, seconds')

for i=1:length(timeSeries)
    monitorTemp(i)=states{i, 1}.T(prod1);
end
monitorTemp=monitorTemp-273.15;
figure,
plot(timeSeries,monitorTemp,'bo','LineWidth',2)
legend('MRST')
ylabel('Temp, C')
xlabel('Time, seconds')

figure,
plotGrid(G,'faceAlpha',0.1,'faceColor','y');
% set(gca,'ZDir','reverse');
% hold on
plotWell(G,W);
axis tight equal;



