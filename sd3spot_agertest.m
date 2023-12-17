function [description, options, state0, model, schedule, plotOptions] = sd3spot_agertest(varargin)
description = 'Small, enhanced geothermal system';
options = struct('ncells',6, 'nstep',60, 'nlayers',9);
options = merge_options(options, varargin{:});
if nargout <= 2, return; end

% Module dependencies
require ad-core ad-props ad-blackoil geothermal compositional upr
gravity reset on


xmaxFS = [1200,600,250];


%working sdf
data.SDF={[ 05, 150;   600,05],...
                [ 05, 150; 600,300],...
                [ 05, 450;  600,300],...
                [ 05, 450;  600,595],...
                [ 600, 05;  1195,150],...
                [ 600, 300;  1195,150],...
                [ 600, 300;  1195,450],...
                [600, 595;  1195,450]};
            

%% 2D fractures and grid processing
physdim= [1,1,1]; % Scaled dimension for proper gridding
xmax=physdim(1:2);
fractureScaled2D= cellfun(@(x) (x./xmaxFS(1:2)), data.SDF, 'UniformOutput', false); % Scaled SDF in 2D


% GRID
dx   = 20./xmaxFS(1); % Approximate cell size
protD = {@(p) 0.00008*ones(size(p,1),1)};
srcRad=2;

G2D = pebiGrid2D(dx, xmax(1:2), ...
    'cellConstraints', fractureScaled2D, ... % Fractures
    'protLayer',true,....
    'protD',protD,...
    'CCRefinement'   , true , ... % Refine fractures
    'CCEps',1/100,...
    'CCFactor' , 0.15  );    % Rel fracture cell size

G2DFD=G2D;
G2DFD.nodes.coords=G2DFD.nodes.coords.*xmaxFS(1:2);



%% Make 3D layered grid
layers = diff(linspace(0, xmaxFS(3), options.nlayers + 1));
G = makeLayeredGrid(G2DFD, layers);
G = computeGeometry(G);

% Convert the fracture to full dim 3D
FractureFullDim=cellfun(@(x) (x.*xmaxFS(1:2)), fractureScaled2D, 'UniformOutput', false);

FractureFulldim3D=[];
for i=1:numel(FractureFullDim)
    Curfracdata=FractureFullDim{1,i};
    FractureFulldim3D=[FractureFulldim3D {[Curfracdata(1,:) 0;Curfracdata(2,:) 0;Curfracdata(2,:) xmaxFS(3);Curfracdata(1,:) xmaxFS(3)]}];
end
G.cells.tag = repmat(G2DFD.cells.tag, options.nlayers, 1);



%% Plot Grids
numSDfrac=numel(data.SDF);


% HR edit: Create Slot Drill Fractures
pts=G.cells.centroids(:,:);
o = Distance(pts,FractureFulldim3D(1:numSDfrac));
nearcell=(o<srcRad);
SdfTagInit=zeros(size(nearcell(:,1)));
for i=1:size(nearcell,2)
    SdfTagInit=SdfTagInit | nearcell(:,i);
end

data.fractures=FractureFulldim3D;
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



for i=1:size(TagAr,2)
    G.cells.SDF=G.cells.SDF | TagAr(:,i);
end

G.cells.tag=G.cells.SDF;




%% Well trajectories
% HR Edit: Co-ordinates for well trajectories
% Vertical production wells (2D coordinates)
fracture_cells = G.cells.tag;
CordProd=[600,05;
          600,595];
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

distances2 = sqrt((x-CordProd(2,1)).^2 + (y-CordProd(2,2)).^2);
[~, distances2] = sort(distances2);


prod1=distances1(1);
prod2=distances2(1);




%Injection well

CordInj=[600,300];
x=G.cells.centroids(fracture_cells,1);
y=G.cells.centroids(fracture_cells,2);

distances1 = sqrt((x-CordInj(1,1)).^2 + (y-CordInj(1,2)).^2);
[~, distances1] = sort(distances1);
Inj=distances1(1); % As it is a vertical well.

% Intersect all the layers
map = find(fracture_cells);
Inj = map(Inj); prod1 = map(prod1); prod2 = map(prod2); 
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
time = 50*year;

%% Set up wells

rw=0.01; % Well radius

W = addWell([], G, rock, Inj, ...
    'type', 'rate', 'val',0.069,'compi', 1, 'name', 'inj','Radius', rw,'Dir','z');
W = addWell(W, G, rock, prod1, ...
    'type', 'bhp', 'val', 25e6, 'compi', 1, 'name', 'prod_1','Radius', rw,'Dir','z');
W = addWell(W, G, rock, prod2, ...
    'type', 'bhp', 'val', 25e6, 'compi', 1, 'name', 'prod_2','Radius', rw,'Dir','z');

W = addThermalWellProps(W, 'T', Tinj);
% Make schedule
schedule = simpleSchedule(rampupTimesteps(time, time/options.nstep,15), 'W', W);
% Plotting
plotOptions = {'View', [65, 30], ...
    'Size', [700, 500], ...
    'PlotBoxAspectRatio', [1, xmaxFS(2)/xmaxFS(1), 15/xmaxFS(1)]};
end
% 
