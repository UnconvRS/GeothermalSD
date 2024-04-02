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

dispersity=1e-9; %larger less dispersion, larger value needed when defined frac_angle
frac_angle=0; %0 means random direction

%% Only SRV

%Connected DFN
%{
fnm=[]
for i=1:NumFracs
    fr=[FracX_Range(1)+(i-1)*Frac_Spacing FracX_Range(1)+(i)*Frac_Spacing];
    Domain=[fr(1),FracY_Range(1),fr(2),FracY_Range(2)];
    set=Field(DFN('dim',2,'n',NumNaturalFracs,'bbx',Domain,...
    'minl',min_frac_length,'mu',avg_frac_length,'maxl',max_frac_length, ...
    'dir',frac_angle,'ddir',-dispersity ...
    ),'Line');
    fnm=[fnm;set];
end

%Resolve intersection split
%[ols,jds,sn,xts,ids,La] = Split(fnm);
%fnm=ols;
%}

%Non-Connected DFN
%{
fnm=[]
for i=1:NumFracs
    fr=[FracX_Range(1)+(i-1)*Frac_Spacing FracX_Range(1)+(i)*Frac_Spacing];
    Domain=[fr(1),40,fr(2),100];
    set=Field(DFN('dim',2,'n',NumNaturalFracs,'bbx',Domain,...
    'minl',min_frac_length,'mu',avg_frac_length,'maxl',max_frac_length, ...
    'asep',0.1),'Line');
    fnm=[fnm;set];
end
%}

%% Whoel Domain
%{
%Domain Size
Domain=[10,10,190,130];
NumNaturalFracs=150;

fnm=Field(DFN('dim',2,'n',NumNaturalFracs,'bbx',Domain,...
    'minl',min_frac_length,'mu',avg_frac_length,'maxl',max_frac_length, ...
    'dir',frac_angle,'ddir',-dispersity ...
    ),'Line');
%}

%%{
%Domain Size
Domain=[10,10,190,130];
NumNaturalFracs=150;

fnm=Field(DFN('dim',2,'n',NumNaturalFracs,'bbx',Domain,...
    'minl',min_frac_length,'mu',avg_frac_length,'maxl',max_frac_length, ...
    'asep',0.1),'Line');
%%}

%Total Length of fracture
d=diff(fnm(1:end,:));
length=sum(sqrt(sum(d.*d,2)));
fprintf('TotalLength=%f\n',length);
Draw('line',fnm);
%restoredefaultpath;