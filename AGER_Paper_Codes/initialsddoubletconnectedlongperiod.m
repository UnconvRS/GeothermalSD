function [description, options, state0, model, schedule, plotOptions] = initialsddoubletconnectedlongperiod(varargin)
description = 'Small, enhanced geothermal system';
options = struct('ncells',6, 'nstep',60, 'nlayers',9);
options = merge_options(options, varargin{:});
if nargout <= 2, return; end

% Module dependencies
require ad-core ad-props ad-blackoil geothermal compositional upr
gravity reset on

Globals;
physdim= [1,1,1]; 
% xmaxFS = [1200,300,250];
xmaxFS = [1200,600,250];
xmax=physdim(1:2);
rng(202095)
tol=1e-5;
%Total Number of Fractures
NumNaturalFracs=0;
Domain=[tol,tol,physdim(1)-tol, ...
    physdim(2)-tol];

%Fractuer Paramters
max_frac_length=.150;
avg_frac_length=.075;
min_frac_length=.025;

set1=Field(DFN('dim',2,'n',NumNaturalFracs,'bbx',Domain,...
    'minl',min_frac_length,'mu',avg_frac_length,'maxl',max_frac_length, ...
    'asep',0.1),'Line');

%working sdf 1
% data.fractures={[ 5, 595;   100,50],...
%                 [ 25, 570;  350,175],...
%                 [ 25, 570;  600,300],...
%                 [ 25, 595;  850,425],...
%                 [ 25, 595;  1100,550],...
%                 [ 100, 50;  1170,5],...
%                 [ 350, 175; 1170,5],...
%                 [600, 300;  1170,25],...
%                 [850, 425;  1170,25],...
%                 [1100, 550;  1195,5]};
            

% % %working sdf 2
% data.fractures={[ 25, 575;   200,100],...
%                 [ 25, 575;  400,200],...
%                 [ 5, 595;  1170,20],...
%                 [ 50, 595;  800,400],...
%                 [ 50, 595;   1000,500],...
%                 [ 200, 100;  1195,5],...
%                 [400, 200; 1170,20],...
%                 [800, 400;  1170,20],...
%                 [1000, 500;  1195,5]};
            
 %working sdf 3
% data.fractures={[ 5, 595;   200,150],...
%                 [ 200, 150;  1195,5],...
%                 [ 5, 595;  600,250],...
%                 [ 600,250;  1195,5],...
%                 [ 5, 595;  900,450],...
%                 [900,450;   1195,5]};    
            
 data.fractures={[ 20, 585;   200,150],...
                [ 200, 150;  1185,20],...
                [ 20, 585;  600,250],...
                [ 600,250;  1185,20],...
                [ 20, 585;  900,450],...
                [900,450;   1185,20]};              
            
 
% %working sdf test
% data.fractures={[ 5, 595;   200,100],...
%                 [ 5, 595;  400,200],...
%                 [ 5, 595;  1195,5],...
%                 [ 5, 595;  800,400],...
%                 [ 5, 595;   1000,500],...
%                 [ 200, 100;  1195,5],...
%                 [400, 200; 1195,5],...
%                 [800, 400;  1195,5],...
%                 [1000, 500;  1195,5]};



data.fractures2=[];
for i=1:numel(data.fractures)
    Curfracdata=data.fractures{1,i};
    data.fractures2=[data.fractures2 {[Curfracdata(1,:) 0;Curfracdata(2,:) 0;Curfracdata(2,:) xmaxFS(3);Curfracdata(1,:) xmaxFS(3)]}];
end
            
FractureScaled= cellfun(@(x) (x-[0,0]), data.fractures, 'UniformOutput', false);
FractureScaled= cellfun(@(x) (x./xmaxFS(1:2)), FractureScaled, 'UniformOutput', false);

FractureScaled2= cellfun(@(x) (x-[0,0,0]), data.fractures2, 'UniformOutput', false);
FractureScaled2= cellfun(@(x) (x./xmaxFS), FractureScaled2, 'UniformOutput', false);
% fracs = FractureScaled';
% for i=1:size(set1,1)
%     fracs=[fracs {[set1(i,1),set1(i,2);set1(i,3),set1(i,4)]}];
% end

fracs2 = [FractureScaled2';set1];


%% GRID
% GRID
dx   = 20./xmaxFS(1); % Approximate cell size
protD = {@(p) 0.00008*ones(size(p,1),1)};
ApertureSD=2;
Aperturenf=0.005;
G2D = pebiGrid2D(dx, xmax(1:2), ...
    'cellConstraints', FractureScaled, ... % Fractures
    'protLayer',true,....
    'protD',protD,...
    'CCRefinement'   , true , ... % Refine fractures
    'CCEps',1/100,...
    'CCFactor' , 0.15  );    % Rel fracture cell size

figure,
plotGrid(G2D,G2D.cells.tag, 'faceColor', [1,1,1]*0.8,'LineWidth',1); axis tight
plotGrid(G2D, 'faceColor', 'none', 'edgeAlpha', 0.1);


G2DFD=G2D;
G2DFD.nodes.coords=G2DFD.nodes.coords.*xmaxFS(1:2);
figure,
plotGrid(G2DFD,G2DFD.cells.tag, 'faceColor', [1,1,1]*0.8,'LineWidth',1); axis tight
plotGrid(G2DFD, 'faceColor', 'none', 'edgeAlpha', 0.1);
%  save('Grid152DFD.mat','G2DFD');

% Make layered grid
layers = diff(linspace(0, xmaxFS(3), options.nlayers + 1));
G = makeLayeredGrid(G2DFD, layers);
G = computeGeometry(G);



% G.cells.tag = repmat(G2DFD.cells.tag, options.nlayers, 1);





%% Plot Grids
FractureFullDim=cellfun(@(x) (x.*xmaxFS), fracs2, 'UniformOutput', false);
numSDfrac=numel(data.fractures);


%Slot Drill Fractures
pts=G.cells.centroids(:,:);
o = Distance(pts,FractureFullDim(1:numSDfrac));
nearcell=(o<ApertureSD);

% Natural fractures

if numel(FractureFullDim)>numSDfrac
    o = Distance(pts,FractureFullDim(numSDfrac+1:end));
    nearcellnf=(o<1.5);
end


% 
G.cells.SDF=zeros(size(nearcell(:,1)));

if numel(FractureFullDim)>numSDfrac
    G.cells.NF=nearcellnf(:,1);
    for i=2:size(nearcellnf,2)
        G.cells.NF=G.cells.NF | nearcellnf(:,i);
    end
else
    G.cells.NF=[];
end

% options.nlayers=ceil(physdim(3)/ds_3(3));

% FractureFullDim2=cellfun(@(x) (x.*xmaxFS), FractureScaled2, 'UniformOutput', false);
% %% Craete slot drill fracture
% data.fractures= cellfun(@(x) (x-[40,40,0]), data.fractures, 'UniformOutput', false);
data.fractures=data.fractures2;
for i=1:numSDfrac
    G.cells.tag1=nearcell(:,i);
    curfracTag=nearcell(:,i);
    fracLen=(sqrt(sum((data.fractures{i}(1,1:2) - data.fractures{i}(2,1:2)) .^ 2)));
    if fracLen> 1200
        FracCentr=[(data.fractures{i}(1,1)+data.fractures{i}(2,1))/2,(data.fractures{i}(1,2)+data.fractures{i}(2,2))/2,data.fractures{i}(1,3)];
        R1= (sqrt(sum((data.fractures{i}(1,1:2) - FracCentr(1:2)) .^ 2)))/2;
        R2= (sqrt(sum((FracCentr(1:2) - data.fractures{i}(2,1:2)) .^ 2)))/2;
        FracCentr1=[(data.fractures{i}(1,1)+FracCentr(1))/2,(data.fractures{i}(1,2)+FracCentr(2))/2,data.fractures{i}(1,3)];
        FracCentr2=[(FracCentr(1)+data.fractures{i}(2,1))/2,(FracCentr(2)+data.fractures{i}(2,2))/2,data.fractures{i}(1,3)];
    
        Dist1=sqrt(sum((G.cells.centroids(curfracTag,:)-FracCentr1).^2,2));
        Dist2=sqrt(sum((G.cells.centroids(curfracTag,:)-FracCentr2).^2,2));
        
        CurcheckDist=Dist1>R1 & Dist2>R2;
    else
        FracCentr=[(data.fractures{i}(1,1)+data.fractures{i}(2,1))/2,(data.fractures{i}(1,2)+data.fractures{i}(2,2))/2,data.fractures{i}(1,3)];
        R= (sqrt(sum((data.fractures{i}(1,1:2) - data.fractures{i}(2,1:2)) .^ 2)))/2;
        Dist=sqrt(sum((G.cells.centroids(curfracTag,:)-FracCentr).^2,2));
        CurcheckDist=Dist>R & G.cells.centroids(curfracTag,3)>(xmaxFS(3)/options.nlayers);
    end
    G.cells.tag1(curfracTag)=~CurcheckDist;
    
    TagAr(:,i)=G.cells.tag1;
end
G.cells.SDF=false(numel(G.cells.tag1),1);
fraclayer1=find(G2DFD.cells.tag);
G.cells.SDF(fraclayer1)=true;
for i=1:size(TagAr,2)
    G.cells.SDF=G.cells.SDF | TagAr(:,i);
end

if sum(G.cells.NF)>1
    G.cells.tag= G.cells.NF | G.cells.SDF;
else
    G.cells.tag=G.cells.SDF;
end

%% Well trajectories
% HR Edit: Co-ordinates for well trajectories
% Vertical production wells (2D coordinates)
fracture_cells = G.cells.tag;
CordProd=[20,585];
%{
        HR Edit:  Finding the well-fracture intersection.
        !!!!  Doesn't work if the well co-ordinates is in the fracture cell
        faces.
        
        Find the eucledian distance from the well co-ordinates to the
        nearest fracture cells.
%}

x=G.cells.centroids(fracture_cells,1);
y=G.cells.centroids(fracture_cells,2);
z=G.cells.centroids(fracture_cells,3);
distances1 = sqrt((x-CordProd(1,1)).^2 + (y-CordProd(1,2)).^2);
[~, distances1] = sort(distances1);



prod1=distances1(1);


%Injection well

CordInj=[1185,20];
x=G.cells.centroids(fracture_cells,1);
y=G.cells.centroids(fracture_cells,2);

distances1 = sqrt((x-CordInj(1,1)).^2 + (y-CordInj(1,2)).^2);
[~, distances1] = sort(distances1);
Inj=distances1(1); % As it is a vertical well.

% Intersect all the layers
map = find(fracture_cells);
Inj = map(Inj); prod1 = map(prod1); 
%%

% Make rock
rock = makeRock(G, 10*nano*darcy, 0.01);
rock.perm(fracture_cells) = 10*darcy; % Fracture permeability
FracCellVolm=G.cells.volumes(fracture_cells);
poreVolm2achieved=4.36;
rock.poro(fracture_cells) = poreVolm2achieved./FracCellVolm; 

% Add thermal properties
Watt = joule/second;
rock = addThermalRockProps(rock, 'lambdaR', 3*Watt/(meter*Kelvin)       , ...
    'rhoR'   , 2700*kilogram/meter^3       , ...
    'CpR'    , 1000*joule/(kilogram*Kelvin));
K0   = 273.15*Kelvin;

% Make fluid
fluid = initSimpleADIFluid(                     ...
    'phases', 'W'                  , ... % Water only
    'mu'    , 0.5*centi*poise      , ... % Viscosity
    'rho'   , 1000*kilogram/meter^3, ... % Reference density
    'c'     , 4.4e-10/Pascal       , ... % Compressibility
    'pRef'  , 1*atm                );    % Ref pressure

% Add thermal properties
fluid = addThermalFluidProps(fluid            , ...
    'Cp'     , 4.2*joule/(gram*Kelvin), ... % Heat capacity
    'lambdaF', 0.6*Watt/(meter*Kelvin), ... % Thermal cond
    'cT'     , 207e-6/Kelvin          , ... % Thermal expansion
    'TRef'   , K0 + 10*Kelvin         );    % Reference temp

% Make model
model = GeothermalModel(G, rock, fluid);
model.extraStateOutput = true;
model.outputFluxes     = true;

% Initial state
% HR Edit: Find the temperature using geothermal gradient
Tres = 490*Kelvin+G.cells.centroids(:,3)*(50/1000); %473  temp at top 50k/km geothermal gradient
Tinj = 293*Kelvin; %60F ambient T
state0   = initResSol(G, 30e6, 1);   %(G, 300*barsa, 1);
state0.T = Tres;
time = 250*year;

%% Set up wells

rw=0.01; % Well radius
% W = addWell([], G, rock, Inj, ...
%     'type', 'bhp', 'val', 40e6,'compi', 1, 'name', 'inj','Radius', rw,'Dir','z');
W = addWell([], G, rock, Inj, ...
    'type', 'rate', 'val',0.069,'compi', 1, 'name', 'inj','Radius', rw,'Dir','z');
% W = addWell(W, G, rock, prod1, ...
%     'type', 'rate', 'val',-0.069,'compi', 1, 'name', 'prod','Radius', rw,'Dir','z');
W = addWell(W, G, rock, prod1, ...
    'type', 'bhp', 'val', 25e6, 'compi', 1, 'name', 'prod_1','Radius', rw,'Dir','z');


W = addThermalWellProps(W, 'T', Tinj);
% Make schedule
schedule = simpleSchedule(rampupTimesteps(time, time/options.nstep,15), 'W', W);
% Plotting
plotOptions = {'View', [65, 30], ...
    'Size', [700, 500], ...
    'PlotBoxAspectRatio', [1, xmaxFS(2)/xmaxFS(1), 15/xmaxFS(1)]};
end
% 
