function [E] = ele_strain_energy(C, Be,struc, U,edofMat1)
nely = size( struc,1); nelx = size( struc, 2); nelz = size( struc,3);
nele = nelx*nely*nelz;
E = zeros(nele,1);
ele_epsilon = Be*U(edofMat1');
sigma = zeros(nele,6); % Stress Components of each element [ xx, yy, zz, xy,yz,xz]

m =0;
for elz = 1:nelz
    for ely = 1:nely
        for elx = 1:nelx
            m = m + 1;
            sigma(m,:) = struc(ely,elx,elz)*( C*ele_epsilon(:,m) )';
            E(m,1) = sigma(m,:)*ele_epsilon(:,m);
        end
    end
end