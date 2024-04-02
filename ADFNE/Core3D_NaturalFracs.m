Globals;
global Tolerance;
rng(1234567890);

%Set up proper tolerance, default 1e-14
DomainSize=[0.5,0.5,0.5,2.5,1.5,1.5];
Tolerance = Tolerance*max(DomainSize);

%Generate the abritary shape domain
plys = GetShapePoly('shape','Cylinder','Radius',0.5,...
                    'Center',[0.5,1,1],'Length',2,'Resolution',20);
Draw('ply',plys,'fa',0.1,'fc',[0, 0.4470, 0.7410]);% Bounding domain polygon mesh
       
%Frac parameters
TotalNumOfFracs=100;
minLength=0.1;
maxLength=0.3;
averageLength=0.2;

fnm = Field(DFN('dim',3,'n',TotalNumOfFracs-1,...
            'minl',minLength,'maxl',maxLength,'mu',averageLength,...
            'bbx',[0.45,0.45,0.45,2.55,1.55,1.55],'polyDomain',plys,...
            'shape','c','q',12),'Poly'); %Large Polygon DFN
       
%Find frac intersections
[xts,ids,La] = Intersect(fnm);

%Plot fracs and its intersections
clf;
Draw('ply',plys,'fa',0.1,'fc',[0, 0.4470, 0.7410]);% Bounding domain polygon mesh
Draw('lin',xts);% 3D DFN Model
hold on;
Draw('ply',fnm);% 3D DFN Model