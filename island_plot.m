figure(10); hold on

for i = 1999-20: 2000-20
    [island_hatch_line_xy,poly_island] = intersect(i, edofMat, nodeCor, theta(i) );
    mapshow( poly_island(1,:), poly_island(2,:),'LineWidth', 2 );
    num_hatch_lines = size(island_hatch_line_xy,1);
    
    for k = 1:2:num_hatch_lines
        mapshow( island_hatch_line_xy(k, [1,3] ), island_hatch_line_xy(k, [2,4] ),'Color','r' );
    end
     
end
