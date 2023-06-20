H = 1000;
Merged_layer = dz/(30/H/H);
Layer_thickness = 30/H; % unit: mm


FileID = fopen( 'Block_PITT.cli', 'w');
fprintf(FileID,'$$HEADERSTART\n');
fprintf(FileID,'$$ASCII\n');
fprintf(FileID,'$$UNITS/1.000000\n');
fprintf(FileID,'$$VERSION/200\n');
fprintf(FileID,'$$HEADEREND\n');
fprintf(FileID,'$$GEOMETRYSTART\n');
fprintf(FileID,'$$LAYER/0.00000\n');

for m = 1:nelz
    
    hatch_line_xy = [];
    for ele = (m-1)*eleIncrement + 1 : m*eleIncrement
        [island_hatch_line_xy,~] = intersect(ele, edofMat, nodeCor, theta(ele) );
        hatch_line_xy = [hatch_line_xy; island_hatch_line_xy];
    end
    
    num_hatch_line = size( hatch_line_xy, 1);
    
    for i = 1:Merged_layer
        fprintf(FileID, strcat('$$LAYER/',num2str( (m-1)*dz*H + i*Layer_thickness,...
            '%.5f'),'\n') );
        fprintf(FileID, strcat('$$HATCHES/1,', num2str(num_hatch_line,'%d'),',') );
        for j = 1 : num_hatch_line - 1
            if mod(j,2) == 1
                fprintf(FileID, strcat( num2str( hatch_line_xy(j,1) ,'%.5f'),',',...
                    num2str( hatch_line_xy(j,2) ,'%.5f'),',',...
                    num2str( hatch_line_xy(j,3) ,'%.5f'),',',...
                    num2str( hatch_line_xy(j,4) ,'%.5f'),','...
                    ) );
            else
                fprintf(FileID, strcat( num2str( hatch_line_xy(j,3) ,'%.5f'),',',...
                    num2str( hatch_line_xy(j,4) ,'%.5f'),',',...
                    num2str( hatch_line_xy(j,1) ,'%.5f'),',',...
                    num2str( hatch_line_xy(j,2) ,'%.5f'),','...
                    ) );                
            end   
        end
        
        if mod( num_hatch_line, 2) == 1
            fprintf(FileID,strcat( num2str( hatch_line_xy(num_hatch_line,1) ,'%.5f'),',',...
                num2str( hatch_line_xy(num_hatch_line,2) ,'%.5f'),',',...
                num2str( hatch_line_xy(num_hatch_line,3) ,'%.5f'),',',...
                num2str( hatch_line_xy(num_hatch_line,3) ,'%.5f'),'\n'...
                ) );
        else
            fprintf(FileID,strcat( num2str( hatch_line_xy(num_hatch_line,3) ,'%.5f'),',',...
                num2str( hatch_line_xy(num_hatch_line,4) ,'%.5f'),',',...
                num2str( hatch_line_xy(num_hatch_line,1) ,'%.5f'),',',...
                num2str( hatch_line_xy(num_hatch_line,2) ,'%.5f'),'\n'...
                ) );            
        end
    end
    
end

fclose(FileID);