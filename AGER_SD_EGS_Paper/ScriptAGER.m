%% Author : Harun Rashid (hrashid.lsu@gmail.com)
% Add necessary MRST modules-
clc;clear;close all;
mrstModule add ad-core ad-props ad-blackoil geothermal compositional upr...
    mrst-gui example-suite
mrstVerbose on

% Set up example
%% Main Text  
%(Uncomment a specific case and comment others to run that specific case)

example = MRSTExample('fig8ato10'); % Fig 8a, 9, 10a-c
% example = MRSTExample('fig8bto10'); % Fig 8b, 9, 10d-f
% example = MRSTExample('fig8cto10'); % Fig 8c, 9, 10g-i

% example = MRSTExample('fig11ato13a'); % Fig 11a, 12, 13a
% example = MRSTExample('fig11bto13b'); % Fig 11b, 12, 13b

% example = MRSTExample('fig14a_16a_17'); % Fig 14a, 16a, 17
% example = MRSTExample('fig14b_16b_17'); % Fig 14b, 16b, 17

% example = MRSTExample('fig15a_17'); % Fig 15a, 17
% example = MRSTExample('fig15b_17'); % Fig 15b, 17

% example = MRSTExample('fig18a_19a_20'); % Fig 18a, 19a, 20
% example = MRSTExample('fig18b_19b_20'); % Fig 18b, 19b, 20

% example = MRSTExample('fig21a_22a_23'); % Fig 21a, 22a, 23
% example = MRSTExample('fig21b_22b_23'); % Fig 21b, 22b, 23

%% Supplementary
% example = MRSTExample('figS1atoS2'); % Fig S1a, S2
% example = MRSTExample('figS3atoS4'); % Fig S3a, S4

% example = MRSTExample('figS5atoS6'); % Fig S5a, S6
% example = MRSTExample('figS5btoS6'); % Fig S5b, S6

%% Plot setup (HR)
 
WellCell=[];
for i=1:numel(example.schedule.control.W)
    WellCell=[WellCell;example.schedule.control.W(i).cells];
end

example.figure();
if isfield(example.model.G.cells,'SDF')
    plotGrid(example.model.G, example.model.G.cells.SDF , 'faceColor', 'y', 'edgeColor', 'none');
elseif isfield(example.model.G.cells,'HF')
    plotGrid(example.model.G, example.model.G.cells.HF , 'faceColor', 'y', 'edgeColor', 'none');
else
    plotGrid(example.model.G, example.model.G.cells.tag , 'faceColor', 'y', 'edgeColor', 'none');
end
if isfield(example.model.G.cells,'NF')
    plotGrid(example.model.G, example.model.G.cells.NF , 'faceColor', 'b', 'edgeColor', 'none');
end
plotGrid(example.model.G,'faceColor', 'none', 'edgeAlpha', 0.1);
plotWell(example.model.G, example.schedule.control(1).W, 'color', 'r', 'fontSize', 18);
example.setAxisProperties(gca);
axis equal;



FracV=sum(example.model.G.cells.volumes(example.model.G.cells.tag));
PoreVolm=sum(example.model.G.cells.volumes(example.model.G.cells.tag).*example.model.rock.poro(example.model.G.cells.tag));
G=example.model.G;
G.cells.wells=false(length(G.cells.tag),1);
G.cells.wells(WellCell)=true;
G.cells.SDFTag(WellCell)=false;

%% Simulate
problem = example.getPackedSimulationProblem();
clearPackedSimulatorOutput(problem, 'prompt', true);
simulatePackedProblem(problem);
[wellSols, states, reports] = getPackedSimulatorOutput(problem);


%% Results/post processing (HR)

figure,
plotCellData(example.model.G,states{27,1}.T,'edgeAlpha',0);
set(gca,'clim',[300,500],'FontName','Arial','fontSize',14,'xtick',[],'ytick',[]);
colormap jet;
axis tight equal;
view(0,90);
colorbar;

figure,
plotCellData(example.model.G,states{51,1}.T,'edgeAlpha',0);
set(gca,'clim',[300,500],'FontName','Arial','fontSize',14,'xtick',[],'ytick',[]);
colormap jet;
axis tight equal;
view(0,90);
colorbar;


figure,
plotCellData(example.model.G,states{end,1}.T,'edgeAlpha',0);
set(gca,'clim',[300,500],'FontName','Arial','fontSize',14,'xtick',[],'ytick',[]);
colormap jet;
axis tight equal;
view(0,90);
colorbar;




figure,
plotToolbar(example.model.G,states);
colormap jet;
axis tight equal;
colorbar;


% Plot well solutions
plotWellSols(wellSols, example.schedule.step.val)


%thermal energy calculation
p   = getWellOutput(wellSols, 'bhp');
T   = getWellOutput(wellSols, 'T');
q   = abs(getWellOutput(wellSols, 'qWs'));

[h, rho] = deal(zeros(size(p)));
if size(h,2)==2
    for i = 1:2
        h(:, i)   = example.model.fluid.hW(p(:,i), T(:,i));
        rho(:, i) = example.model.fluid.rhoW(p(:,i), T(:,i));
    end
elseif size(h,2)==3
    for i = 1:3
        h(:, i)   = example.model.fluid.hW(p(:,i), T(:,i));
        rho(:, i) = example.model.fluid.rhoW(p(:,i), T(:,i));
    end
end

pref=1*atm;
Tref=273.15;
href= example.model.fluid.hW(pref, Tref);


qH  = abs(q.*rho.*(h-href));

if size(qH,2)==2
    eff = qH(:,2)./qH(:,1);
elseif size(qH,2)==3
    eff = (qH(:,2)+qH(:,3))./qH(:,1);
end

nsteps=size(example.schedule.step.val,1);
time = cumsum(example.schedule.step.val);
if size(qH,2)==2
    qwhCume=cumtrapz(time,qH(:,2));
elseif size(qH,2)==3
    qwhCume=cumtrapz(time,qH(:,2)+qH(:,3));
end
figure,
plot(time/year, qwhCume, 'color', 'k', 'linew', 2);
% axis([[time(5), time(end)]/year, min(eff(5:end))*0.95, max(eff(5:end))*1.05]);
set(gca, 'Box', true, 'FontSize', 13);
xlabel('Time (years)')
title('Cumulative produced energy')






% Calculation of recovery fraction (HR)
G=example.model.G;
tol=1e-3;
Vtotal=sum(G.cells.volumes);
TrO=298;
TrI=mean(example.state0.T);
for i=1:nsteps
    TrA=sum(states{i}.T.*G.cells.volumes)./sum(G.cells.volumes);
    activeCell= (states{i}.T+tol<example.state0.T);
    Vactive=sum(G.cells.volumes(activeCell));
    Fr(i)=(Vactive/Vtotal)*((TrI-TrA)/(TrI-TrO));
end
figure,
plot(time/year,Fr);

%Calculation of Average Thermal Energy
AvgTempRes=sum(states{75, 1}.T.*G.cells.volumes)/(sum(G.cells.volumes));