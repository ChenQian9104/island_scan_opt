figure(4)
m = 0;
for k = 1:nelz
    for j = 1:nely
        for i = 1:nelx
            m = m + 1;
            if j >= 0
                dofs = edofMat1(m, 3:3:24)';
                dis = 1000*U(dofs,:);
                vert = 1000*[ nodeCor( edofMat(m,1),: ); nodeCor( edofMat(m,2),: );... 
                    nodeCor(edofMat(m,3),:); nodeCor(edofMat(m,4),:); ...
                    nodeCor(edofMat(m,5),:); nodeCor(edofMat(m,6),:); ...
                    nodeCor(edofMat(m,7),:); nodeCor(edofMat(m,8),:)];
                face = [ 1 2 3 4; 3 2 6 7; 7 6 5 8; 8 5 1 4; 3 4 8 7; 1 2 6 5];
                patch('Faces',face,'Vertices',vert,'CData',dis,'FaceColor','flat');
                colormap('jet'); 
                axis equal; axis tight;
                view([-45 30]); hold on;
            end
        end
    end
end
