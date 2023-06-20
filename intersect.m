function [hatch_line_xy, poly_island] = intersect( num_island, edofMat, nodeCor, theta)

H = 1000;
hatch = 100/H; % hatching space
poly_island = [];
poly_island_vertice = edofMat(num_island, 1:4);

RotGlob2Loc = zeros(2,2);
RotGlob2Loc(1,1) = cos(theta); RotGlob2Loc(1,2) = sin(theta);
RotGlob2Loc(2,1) = -sin(theta); RotGlob2Loc(2,2) = cos(theta);

RotLoc2Glob = zeros(2,2);
RotLoc2Glob(1,1) = cos(theta); RotLoc2Glob(1,2) = -sin(theta);
RotLoc2Glob(2,1) = sin(theta); RotLoc2Glob(2,2) = cos(theta);



for i = 1:4
    poly_island = [ poly_island, H*[ nodeCor( edofMat( num_island, i), 1);...
        nodeCor( edofMat( num_island, i), 2)] ];
end
poly_island = [ poly_island, H*[ nodeCor( edofMat( num_island, 1), 1);...
    nodeCor( edofMat( num_island, 1), 2)] ];

poly_island_loc = RotGlob2Loc*poly_island;

xMin = min(poly_island_loc(1,:)); xMax = max( poly_island_loc(1,:) );
yMin = min(poly_island_loc(2,:)); yMax = max( poly_island_loc(2,:) );

x1 = [];y1 = [];
for ix = 1:ceil( (yMax-yMin)/hatch ) 
    x1 = [x1; nan;  xMin - 1; xMax + 1];
    y1 = [y1; nan; yMin + hatch/2 + (ix-1)*hatch; yMin + hatch/2 + (ix-1)*hatch];
end

[xi,yi,~] = polyxpoly(poly_island_loc(1,:)',poly_island_loc(2,:)',x1,y1);

hatch_line =  RotLoc2Glob*[xi';yi'];  % glob coordinate
%mapshow( poly_island(1,:),  poly_island(2,:), 'DisplayType','polygon');
%hold on;
n_hatch_line = size(hatch_line,2)/2;
hatch_line_xy = [];
for i = 1:n_hatch_line
    %mapshow( hatch_line(1,2*i-1:2*i), hatch_line(2,2*i-1:2*i));
    % X1 Y1 X2 Y2
    hatch_line_xy = [hatch_line_xy; hatch_line(1, 2*i-1),hatch_line(2, 2*i-1), hatch_line(1, 2*i),hatch_line(2, 2*i)];
end

end

