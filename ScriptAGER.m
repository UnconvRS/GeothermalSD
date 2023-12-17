% SD Enhanced geothermal system (EGS)
%% Add necessary MRST modules-
clc;clear;close all;
mrstModule add ad-core ad-props ad-blackoil geothermal compositional upr...
    mrst-gui example-suite
mrstVerbose on

%% Set up example
% Initial Cases
example = MRSTExample('sd3spot_agertest');


WellCell=[];
for i=1:numel(example.schedule.control.W)
    WellCell=[WellCell;example.schedule.control.W(i).cells];
end

% Plot Grid
example.figure();
if isfield(example.model.G.cells,'SDF')
    plotGrid(example.model.G, example.model.G.cells.SDF , 'faceColor', 'y', 'edgeColor', 'none');
elseif isfield(example.model.G.cells,'HF')
    plotGrid(example.model.G, example.model.G.cells.HF , 'faceColor', 'y', 'edgeColor', 'none');
% else
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

Cell2Plot=example.model.G.cells.centroids(:,3)>150;
for i=1:length(states)
    states_mod{i, 1}.T=states{i, 1}.T(Cell2Plot);
end


figure,
plotCellData(example.model.G,states_mod{25,1}.T,Cell2Plot,'edgeAlpha',0);
hold on
plotCellData(example.model.G,states{25,1}.T(G.cells.tag),G.cells.tag,'edgeAlpha',0);
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
set(gca, 'Box', true, 'FontSize', 13);
xlabel('Time (years)')
title('Cumulative produced energy')



% Calculation of recovery fraction
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
