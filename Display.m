%% DISPLAY THE L_BRACKET
figure(6); hold on;
for i =1:nele
    vert = 1000*[ nodeCor(edofMat(i,1),:); nodeCor(edofMat(i,2),:);...
        nodeCor(edofMat(i,3),:); nodeCor(edofMat(i,4),:); ...
        nodeCor(edofMat(i,5),:); nodeCor(edofMat(i,6),:); ...
        nodeCor(edofMat(i,7),:); nodeCor(edofMat(i,8),:)];
    face = [ 1 2 3 4; 3 2 6 7; 7 6 5 8; 8 5 1 4; 3 4 8 7; 1 2 6 5];
    %s=patch('Faces', face, 'Vertices',vert,'FaceColor',[30*VM_stress(i),0,0] );
    patch('Faces',face,'Vertices',vert,'CData',[ 1, 1,1,...
        1,1,1],'FaceColor','flat');
   
    axis equal; axis tight;box on;view([-45 45]);
end
