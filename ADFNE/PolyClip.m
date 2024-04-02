function [y,b] = PolyClip(x,plys)
% Clip
% clips "x", polygons [3D] by a closed polygon
%
% Usage...:
% [y,b] = PolyClip(x,plys);
% 
% Input...: x         (n,cell),lines|polygons
%           plys      (n,cell), closed polygons
% Output..: y         (k,cell), clipped
%           b         (n), boolean,record of those fully clipped out
%
% Examples:
%{
[y,b] = Clip(rand(10,4),[0,0,1,1]); % clipping 2d lines
[y,b] = Clip(plys,[0,0,0,1,1,1]); % clipping 3d polygons
%}
%
% Author: Bin Wang
% Updated.: 2019-08-08
if iscell(x)                                                                    % 3D polygons
    nClip = size(plys,1);
    n = size(x,1);
    y = cell(n,1);
    k = 0;                                                                      % counter
    b = false(n,1);                                                             % mask for all lines|polygons
    for i = 1:n
        ply = x{i};
        for ci=1:nClip
            ClipPlane= Plane(plys{ci},true);
            ply = CExt.ClipPolyHP(ply,ClipPlane);
            %if isempty(ply); continue; end
            %ply_last = ply;
            %Draw('ply',plys{ci});
            %Pts_cntr=Center(plys{ci});
            %Draw('ply',ply_last,'fc',[1,0.5,0.1]);
            %norm= CExt.planeNormal(Plane(plys{ci},true));
            %quiver3(Pts_cntr(1), Pts_cntr(2), Pts_cntr(3), norm(1), norm(2), norm(3));
        end
        if isempty(ply); continue; end
        k = k+1;
        y{k} = ply;                                                             % stores clipped item
        b(i) = true;
    end
    y = y(1:k);
end
end
