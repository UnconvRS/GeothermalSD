function [description, options, state0, model, schedule, plotOptions] = fig11bto13b(varargin)
%% Author : Harun Rashid (hrashid.lsu@gmail.com)

description = 'Small, enhanced geothermal system';
options = struct('ncells',6, 'nstep',60, 'nlayers',9);
options = merge_options(options, varargin{:});
if nargout <= 2, return; end

% Module dependencies
require ad-core ad-props ad-blackoil geothermal compositional upr
gravity reset on


physdim= [1,1,1]; 
% xmaxFS = [1200,300,250];
xmaxFS = [1200,600,250];
xmax=physdim(1:2);

rng(202095)
tol=1e-3;
Globals;
set1 = Field(DFN('dim',3,'n',0,'dir',15,'ddir',-100,'minl',0.01,...
    'mu',0.05,'maxl',0.1,'bbx',[tol,tol,tol,physdim(1)-tol, ...
    physdim(2)-tol,physdim(3)-tol],'dip',45,'ddip',-100,...
    'shape','l','q',4),'Poly');

startY=0;
EndY=600;


fractureSpacing=130;
x1=80; %Start point
x2=x1+fractureSpacing;
x3=x2+fractureSpacing;
x4=x3+fractureSpacing;
x5=x4+fractureSpacing;
x6=x5+fractureSpacing;
x7=x6+fractureSpacing;
x8=x7+fractureSpacing;
x9=x8+fractureSpacing;

data.fractures={[x1,startY;x1,EndY],...
    [x2,startY;x2,EndY],...
    [x3,startY;x3,EndY],...
    [x4,startY;x4,EndY],...
    [x5,startY;x5,EndY],...
    [x6,startY;x6,EndY],...
    [x7,startY;x7,EndY],...
    [x8,startY;x8,EndY],...
    [x9,startY;x9,EndY]}; 

data.fractures2=[];
for i=1:numel(data.fractures)
    Curfracdata=data.fractures{1,i};
    data.fractures2=[data.fractures2 {[Curfracdata(1,:) 0;Curfracdata(2,:) 0;Curfracdata(2,:) xmaxFS(3);Curfracdata(1,:) xmaxFS(3)]}];
end
%             
% 
FractureScaled= cellfun(@(x) (x-[0,0]), data.fractures, 'UniformOutput', false);
FractureScaled= cellfun(@(x) (x./xmaxFS(1:2)), FractureScaled, 'UniformOutput', false);

FractureScaled2= cellfun(@(x) (x-[0,0,0]), data.fractures2, 'UniformOutput', false);
FractureScaled2= cellfun(@(x) (x./xmaxFS), FractureScaled2, 'UniformOutput', false);
fracs = [FractureScaled';set1];
fracs2 = [FractureScaled2';set1];


%% Plot Grids
FractureFullDim=cellfun(@(x) (x.*xmaxFS), fracs2, 'UniformOutput', false);
numHfrac=numel(data.fractures);


NumStages=numHfrac;
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





%Slot Drill Fractures
pts=G.cells.centroids(:,:);
o = Distance(pts,FractureFullDim(1:numHfrac));
nearcell=(o<ApertureSD);

% Natural fractures

if numel(FractureFullDim)>numHfrac
    o = Distance(pts,FractureFullDim(numHfrac+1:end));
    nearcellnf=(o<1.5);
end


% 

G.cells.HF=nearcell(:,1);
for i=2:size(nearcell,2)
    G.cells.HF=G.cells.HF | nearcell(:,i);
end
if numel(FractureFullDim)>numHfrac
    G.cells.NF=nearcellnf(:,1);
    for i=2:size(nearcellnf,2)
        G.cells.NF=G.cells.NF | nearcellnf(:,i);
    end
else
    G.cells.NF=[];
end

% options.nlayers=ceil(physdim(3)/ds_3(3));



if sum(G.cells.NF)>1
    G.cells.tag= G.cells.NF | G.cells.HF;
else
    G.cells.tag=G.cells.HF;
end

 fracture_cells = G.cells.tag;
    
 %% Well trajectories
 % HR Edit: Co-ordinates for well trajectories
 % 
    perc=0; %0%    
%     perc=31.25; %25%
%     perc=62.5; %50%
%     perc=93.75; %75%
    wyloc=(startY+EndY)/2;
    CordInj=[x1,wyloc,250-perc;
        x2,wyloc,250-perc;
        x3,wyloc,250-perc;
        x4,wyloc,250-perc;
        x5,wyloc,250-perc;
        x6,wyloc,250-perc;
        x7,wyloc,250-perc;
        x8,wyloc,250-perc;
        x9,wyloc,250-perc]; 
    
    inj=[];
    for i=1:size(CordInj,1)
        distances1=sqrt(sum((G.cells.centroids(fracture_cells,:)-CordInj(i,:)).^2,2));
        [~, distances1] = sort(distances1,'ascend');
        inj(i)=distances1(1);
    end

    
    %Production well
    ProdInj=[x1,wyloc,0+perc;
        x2,wyloc,0+perc;
        x3,wyloc,0+perc;
        x4,wyloc,0+perc;
        x5,wyloc,0+perc;
        x6,wyloc,0+perc;
        x7,wyloc,0+perc;
        x8,wyloc,0+perc;
        x9,wyloc,0+perc];
    
       % To Do:   Find the CordInj and ProdInj automatically
       prod=[];
       for i=1:size(ProdInj,1)
           distances1=sqrt(sum((G.cells.centroids(fracture_cells,:)-ProdInj(i,:)).^2,2));
           [~, distances1] = sort(distances1,'ascend');
           prod(i)=distances1(1);
       end
       
       map = find(fracture_cells);
       inj = map(inj); prod = map(prod);
       
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
    Tres = 490*Kelvin+G.cells.centroids(:,3)*(50/1000); %490  temp at top 50k/km geothermal gradient
    Tinj = 293*Kelvin; %60F ambient T
    state0   = initResSol(G, 30e6, 1);   %(G, 300*barsa, 1);  
    state0.T = Tres;
    time = 50*year;

    %% Make wells
    rw=0.01; % Well radius
    W = addWell([], G, rock, inj, ...
        'type', 'rate', 'val',0.069, 'compi', 1, 'name', 'inj','Radius', rw,'Dir','x');
    W = addWell(W, G, rock, prod, ...
        'type', 'bhp', 'val', 25e6, 'compi', 1, 'name', 'prod_1','Radius', rw,'Dir','x');
    W = addThermalWellProps(W, 'T', Tinj);
    % Make schedule
    schedule = simpleSchedule(rampupTimesteps(time, time/options.nstep,15), 'W', W);
    % Plotting
    plotOptions = {'View', [65, 30], ...
                   'Size', [700, 500], ...
                   'PlotBoxAspectRatio', [1, xmaxFS(2)/xmaxFS(1), 15/xmaxFS(1)]};
end