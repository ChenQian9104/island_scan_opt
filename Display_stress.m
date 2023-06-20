vm = zeros(nele,1);    % Von-Mises stress of each element
sigma = zeros(nele,6); % Stress Components of each element [ xx, yy, zz, xy,yz,xz]
ele_epsilon = Be*U(edofMat1');
m = 0;
element_to_show = [];
for k = 1:nelz
    element_to_show =[ element_to_show,...
        (k-1)*eleIncrement + 1: 1: (k-1)*eleIncrement  + nelx];
end



for elz = 1:nelz
    for ely = 1:nely
        for elx = 1:nelx
            m = m + 1;
            sigma(m,:) = ( C*( ele_epsilon(:,m) - epsilon(:,m ))';
            vm(m,1) = ( sigma(m,1) - sigma(m,2) )^2 + ( sigma(m,2) - sigma(m,3) )^2 + ...
                ( sigma(m,3) - sigma(m,1) )^2 + 6*( sigma(m,4)^2 + sigma(m,5)^2 + sigma(m,6)^2 ); 
            vm(m,1) = sqrt( vm(m,1)/2 );
        end
    end
end

for i = 1: size(element_to_show,2)
    m = element_to_show(i);
    vert = [ nodeCor( edofMat(m,1),: ); nodeCor( edofMat(m,2),: );...  
        nodeCor(edofMat(m,3),:); nodeCor(edofMat(m,4),:); ...
        nodeCor(edofMat(m,5),:); nodeCor(edofMat(m,6),:); ...
        nodeCor(edofMat(m,7),:); nodeCor(edofMat(m,8),:)];
    face = [ 1 2 3 4; 3 2 6 7; 7 6 5 8; 8 5 1 4; 3 4 8 7; 1 2 6 5];
    patch('Faces',face,'Vertices',vert,'CData',[ vm(m), vm(m),vm(m),...
        vm(m),vm(m),vm(m)],'FaceColor','flat'); 
    colormap('jet'); 
    axis equal; axis tight;
    box on;view([-45 45]); 
end




for i = 1: nele
    m = i
    vert = [ nodeCor( edofMat(m,1),: ); nodeCor( edofMat(m,2),: );...  
        nodeCor(edofMat(m,3),:); nodeCor(edofMat(m,4),:); ...
        nodeCor(edofMat(m,5),:); nodeCor(edofMat(m,6),:); ...
        nodeCor(edofMat(m,7),:); nodeCor(edofMat(m,8),:)];
    face = [ 1 2 3 4; 3 2 6 7; 7 6 5 8; 8 5 1 4; 3 4 8 7; 1 2 6 5];
    patch('Faces',face,'Vertices',vert,'CData',[ sigma(m,1), sigma(m,1),sigma(m,1),...
       sigma(m,1),sigma(m,1),sigma(m,1)],'FaceColor','flat'); 
    colormap('jet'); 
    axis equal; axis tight;
    box on;view([-45 45]); 
end


