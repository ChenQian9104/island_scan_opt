% element_to_show = [];
% for k = 1:nelz
%     element_to_show =[ element_to_show,...
%         (k-1)*eleIncrement + nelx*nely/2 + 1: 1: (k-1)*eleIncrement + nelx*nely/2 + nelx];
% end
%element_to_show = setdiff(element_to_show, element_cutting);
%% calculate strain energy of each element
Energy = zeros(nele,1);
ele_epsilon = Be*U(edofMat1');
sigma = zeros(nele,6); % Stress Components of each element [ xx, yy, zz, xy,yz,xz]
m =0;
for elz = 1:nelz
    for ely = 1:nely
        for elx = 1:nelx
            m = m + 1;
            sigma(m,:) = ( C*( ele_epsilon(:,m) - ele_epsilon_act(:,m ) - epsilon(:,m) ) )';
            Energy(m,1) = 0.5*sigma(m,:)*ele_epsilon(:,m);

        end
    end
end



dx = 1; dy = 1; dz = 1;
m = 0;
for j = 1:nely + 1
    for i = 1:nelx + 1
        m = m + 1;
        nodeCor(m,:) = [(i-1)*dx, (j-1)*dy,0];
    end
end

for k = 1:nelz
    for i = 1:nodeIncrement
        nodeCor( i + k*nodeIncrement,:) = nodeCor( i + (k-1)*nodeIncrement,:) + [ 0 0 dz];
    end
end



figure(1)
for i = 1: nele
    m = i
    vert = [ nodeCor( edofMat(m,1),: ); nodeCor( edofMat(m,2),: );...  
        nodeCor(edofMat(m,3),:); nodeCor(edofMat(m,4),:); ...
        nodeCor(edofMat(m,5),:); nodeCor(edofMat(m,6),:); ...
        nodeCor(edofMat(m,7),:); nodeCor(edofMat(m,8),:)];
    face = [ 1 2 3 4; 3 2 6 7; 7 6 5 8; 8 5 1 4; 3 4 8 7; 1 2 6 5];
    patch('Faces',face,'Vertices',vert,'CData',[ Energy(m), Energy(m),Energy(m),...
        Energy(m),Energy(m),Energy(m)],'FaceColor','flat'); 
    colormap('jet'); 
    axis equal; axis tight;
    box on;view([-45 45]); 
end


