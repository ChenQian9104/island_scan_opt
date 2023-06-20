function dv = sensitivity(Be, C, depsilon, U, edofMat, nelx, nely,nelz, dx, dy, dz)
nele = nelx*nely*nelz;
dv = zeros(nele,1);
for i = 1:nele
    fth_dtheta = thermal_load(Be,C,depsilon(:,i),dx,dy,dz);
    dv(i) = fth_dtheta'*U( edofMat(i,:)' );
end

end 