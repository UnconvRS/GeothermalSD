function [plys] = GetShapePoly(varargin)
% GetShapePoly
% generates standard polygon shapes, 3D
%
% Usage...:
% out = DFN(varargin);
%
% Input...: varargin  any
% Output..: out       struct,{Poly,Orig}
%
% Examples:
%{
plys = GetShapePoly('shape','Cylinder','Radius',0.5,...
                    'Center',[1,0.5,0.5],'Length',3,'Resolution',10);      %3D Cylinder
%}
%
% Author: Bin Wang
% Updated.: 2019-08-08
opt = Option(varargin,'shape','Cylinder','Radius',100,'Center',[0,0,0],...
    'Length',0.5,'Resolution',6 ...              % default arguments
    );  % domain added

switch opt.shape
    case {'Cylinder','cylinder'}                                           % axis on X, plane on YZ
        npoly = 2+opt.Resolution;
        plys  = cell(npoly,1);                                             % initializes polygons
        %Start Plane Circle
        angles=linspace(0,2*pi,opt.Resolution+1);
        CirclePts=[cos(angles)*opt.Radius; sin(angles)*opt.Radius]';
        
        plys{1}=[zeros(opt.Resolution+1,1) flip(CirclePts)];
        plys{2}=[ones(opt.Resolution+1,1).*opt.Length CirclePts];
        for i=1:opt.Resolution
            ply=[0.0,CirclePts(i,1),CirclePts(i,2);
                 0.0,CirclePts(i+1,1),CirclePts(i+1,2);
                 opt.Length,CirclePts(i+1,1),CirclePts(i+1,2);
                 opt.Length,CirclePts(i,1),CirclePts(i,2);
                 0.0,CirclePts(i,1),CirclePts(i,2)];
            plys{2+i}=ply;
        end
        plys = Translate(plys,opt.Center);                                 % 3d polygon
end
%clf; Draw('ply',plys,'axes',true);
end