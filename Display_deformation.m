%% Plot displacement along x direction
figure(1)
for i = 1: nele
    m = i
    
    dofs = edofMat1(m, 1:3:24)';
    dis = U(dofs,:);
    
    vert = [ nodeCor( edofMat(m,1),: ); nodeCor( edofMat(m,2),: );...  
        nodeCor(edofMat(m,3),:); nodeCor(edofMat(m,4),:); ...
        nodeCor(edofMat(m,5),:); nodeCor(edofMat(m,6),:); ...
        nodeCor(edofMat(m,7),:); nodeCor(edofMat(m,8),:)];
    face = [ 1 2 3 4; 3 2 6 7; 7 6 5 8; 8 5 1 4; 3 4 8 7; 1 2 6 5];
    patch('Faces',face,'Vertices',vert,'CData',dis,'FaceColor','flat'); 
    colormap('jet'); 
    axis equal; axis tight;
    box on;view([-45 30]); 
end


figure(2)
for i = 1: nele
    m = i
    
    dofs = edofMat1(m, 2:3:24)';
    dis = U(dofs,:);
    
    vert = [ nodeCor( edofMat(m,1),: ); nodeCor( edofMat(m,2),: );...  
        nodeCor(edofMat(m,3),:); nodeCor(edofMat(m,4),:); ...
        nodeCor(edofMat(m,5),:); nodeCor(edofMat(m,6),:); ...
        nodeCor(edofMat(m,7),:); nodeCor(edofMat(m,8),:)];
    face = [ 1 2 3 4; 3 2 6 7; 7 6 5 8; 8 5 1 4; 3 4 8 7; 1 2 6 5];
    patch('Faces',face,'Vertices',vert,'CData',dis,'FaceColor','flat'); 
    colormap('jet'); 
    axis equal; axis tight;
    box on;%view([-45 45]); 
    view([-45 30]);
end


figure(4)
for i = 1:nele
    m = i
    
    
    dofs = edofMat1(m, 3:3:24)';
    
    dis = U(dofs,:);
    vert = [ nodeCor( edofMat(m,1),: ); nodeCor( edofMat(m,2),: );...  
        nodeCor(edofMat(m,3),:); nodeCor(edofMat(m,4),:); ...
        nodeCor(edofMat(m,5),:); nodeCor(edofMat(m,6),:); ...
        nodeCor(edofMat(m,7),:); nodeCor(edofMat(m,8),:)];
    face = [ 1 2 3 4; 3 2 6 7; 7 6 5 8; 8 5 1 4; 3 4 8 7; 1 2 6 5];
    patch('Faces',face,'Vertices',vert,'CData',dis,'FaceColor','flat');
   
    colormap('jet'); 
    axis equal; axis tight;
    box on;%view([-45 45]); 
    view([-45 30]); hold on;
end



%% plot epislon_x
figure(4)
ele_epsilon = Be*U(edofMat1');
for i = 1: nele
    m = i
    
    vert = [ nodeCor( edofMat(m,1),: ); nodeCor( edofMat(m,2),: );...  
        nodeCor(edofMat(m,3),:); nodeCor(edofMat(m,4),:); ...
        nodeCor(edofMat(m,5),:); nodeCor(edofMat(m,6),:); ...
        nodeCor(edofMat(m,7),:); nodeCor(edofMat(m,8),:)];
    face = [ 1 2 3 4; 3 2 6 7; 7 6 5 8; 8 5 1 4; 3 4 8 7; 1 2 6 5];
    patch('Faces',face,'Vertices',vert,'CData',[ele_epsilon(3,m),ele_epsilon(3,m),ele_epsilon(3,m),...
        ele_epsilon(3,m),ele_epsilon(3,m),ele_epsilon(3,m),ele_epsilon(3,m),ele_epsilon(3,m)],'FaceColor','flat'); 
    colormap('jet'); 
    axis equal; axis tight;
    box on;%view([-45 45]); 
    view([-45 30]);
end

%edge_node = [ nodeIncrement*nelz + nelx + 1: nelx + 1: nnode]';
edge_node = [ nodeIncrement*nelz + (nelx + 1)*nely/2 + 1: ...
    nodeIncrement*nelz + (nelx + 1)*nely/2 + nelx + 1]';
obj_U = U(3*edge_node,:);
plot(1:nelx + 1, obj_U);
