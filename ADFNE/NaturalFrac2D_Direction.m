Globals;
%rng(1234567890); %This make fracture is more predictable

%Fracture parameters
NumFracs=5;
Frac_Spacing=25;
Frac_halfLength=30;

FracX_Range=[50-Frac_Spacing/2 150+Frac_Spacing/2];
FracY_Range=[40 100];

%Total Number of Fractures
NumNaturalFracs=25;

%Fractuer Paramters
max_frac_length=25;
avg_frac_length=10;
min_frac_length=5;

dispersity=1.5e1; %larger less dispersion, larger value needed when defined frac_angle
frac_angle=45; %0 means random direction

%% Whoel Domain
%%{
%Domain Size
Domain=[10,10,190,130];
NumNaturalFracs=150;

fnm=Field(DFN_bin('dim',2,'n',NumNaturalFracs,'bbx',Domain,...
    'minl',min_frac_length,'mu',avg_frac_length,'maxl',max_frac_length, ...
    'dir',frac_angle,'ddir',-dispersity ...
    ),'Line');
%%}

%Total Length of fracture
d=diff(fnm(1:end,:));
length=sum(sqrt(sum(d.*d,2)));
fprintf('TotalLength=%f\n',length);
Draw('line',fnm);
%restoredefaultpath;