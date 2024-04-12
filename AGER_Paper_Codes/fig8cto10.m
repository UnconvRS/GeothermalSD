function [description, options, state0, model, schedule, plotOptions] = fig8cto10(varargin)
 description = 'Small, enhanced geothermal system';
options = struct('ncells',6, 'nstep',60, 'nlayers',9);
options = merge_options(options, varargin{:});
if nargout <= 2, return; end

% Module dependencies
require ad-core ad-props ad-blackoil geothermal compositional upr
gravity reset on

xmaxFS = [1200,600,250];





startY=55;
EndY=545;


fractureSpacing=200;
x1=200; %Start point
x2=x1+fractureSpacing;
x3=x2+fractureSpacing;
x4=x3+fractureSpacing;
x5=x4+fractureSpacing;

data.SDF={[x1,startY;x1,EndY],...
    [x2,startY;x2,EndY],...
    [x3,startY;x3,EndY],...
    [x4,startY;x4,EndY],...
    [x5,startY;x5,EndY]}; 

%% 2D fractures and grid processing
physdim= [1,1,1]; % Scaled dimension for proper gridding
xmax=physdim(1:2);
sdfScaled2D= cellfun(@(x) (x./xmaxFS(1:2)), data.SDF, 'UniformOutput', false); % Scaled SDF in 2D

% Include scaled 2D natural fractures
Globals;
rng(20200629)
%Domain Size
Domain=[10./xmaxFS(1),10./xmaxFS(2),xmax(1),xmax(2)];
%Total Number of Fractures

% NumNaturalFracs=377;
NumNaturalFracs=0;
% NumNaturalFracs=160;

%Fractuer Paramters
max_frac_length=0.2;
avg_frac_length=0.07;
min_frac_length=0.001;

set1=Field(DFN('dim',2,'n',NumNaturalFracs,'bbx',Domain,...
    'minl',min_frac_length,'mu',avg_frac_length,'maxl',max_frac_length, ...
    'asep',0.1),'Line');

% Add NF to SDF
fractureScaled2D=sdfScaled2D;
for i=1:size(set1,1)
    fractureScaled2D=[fractureScaled2D {[set1(i,1),set1(i,2);set1(i,3),set1(i,4)]}];
end


%% GRID
dx   = 20./xmaxFS(1); % Approximate cell size
protD = {@(p) 0.00008*ones(size(p,1),1)};
ApertureSD=2;
Aperturenf=0.0005;
G2D = pebiGrid2D(dx, xmax(1:2), ...
    'cellConstraints', fractureScaled2D, ... % Fractures
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
% FractureFullDim=cellfun(@(x) (x.*xmaxFS), fracs2, 'UniformOutput', false);
numSDfrac=numel(data.SDF);


%Slot Drill Fractures
pts=G.cells.centroids(:,:);
o = Distance(pts,FractureFulldim3D(1:numSDfrac));
nearcell=(o<ApertureSD);
SdfTagInit=zeros(size(nearcell(:,1)));
for i=1:size(nearcell,2)
    SdfTagInit=SdfTagInit | nearcell(:,i);
end


% Natural fractures
TempTag=G.cells.tag;
[rmvID,ID]=find(SdfTagInit);
TempTag(rmvID)=0;
NFtagInit=TempTag;


G.cells.NF=NFtagInit;
G.cells.SDF=SdfTagInit;



for i=1:length(G.cells.NF)
    if G.cells.centroids(i,1)>600 && G.cells.centroids(i,2)>300 
        if G.cells.centroids(i,3)>250 || G.cells.centroids(i,3)<200
            G.cells.NF(i)=0;
        end
    end
    if G.cells.centroids(i,1)>600 && G.cells.centroids(i,2)<300
        if G.cells.centroids(i,3)>200 || G.cells.centroids(i,3)<100
            G.cells.NF(i)=0;
        end
    end
    
    if G.cells.centroids(i,1)<=600 && G.cells.centroids(i,2)>300 
        if G.cells.centroids(i,3)>200 || G.cells.centroids(i,3)<150
            G.cells.NF(i)=0;
        end
    end
    
    if G.cells.centroids(i,1)<=600 && G.cells.centroids(i,2)<300
        if G.cells.centroids(i,3)>75 || G.cells.centroids(i,3)<50
            G.cells.NF(i)=0;
        end
    end
end



if sum(G.cells.NF)>1
    G.cells.tag= G.cells.NF | G.cells.SDF;
else
    G.cells.tag=G.cells.SDF;
end
fracture_cells=G.cells.tag;
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
        x5,wyloc,250-perc]; 
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
        x5,wyloc,0+perc];
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
        % if porosity constant
%     rock.poro(fracture_cells)=0.5;
%     G.cells.volumes(fracture_cells) = poreVolm2achieved./rock.poro(fracture_cells);
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
%     Tres=(220+273.15)*ones(G.cells.num,1);
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