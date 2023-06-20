function [epsilon,depsilon] = Inherent_strain( theta )     % Output the inherent strains of each element
nely = size(theta,1); nelx = size(theta,2); nelz = size(theta,3); 
epsilon = zeros(6,nelx*nely*nelz); depsilon = zeros(6,nelx*nely*nelz);

%strain_x = -0.0145; strain_y = -0.0065; strain_z = 0.012;
strain_x = -0.02; strain_y = -0.01; strain_z = 0.015;
strain_xy = sqrt( strain_x^2 + strain_y^2 ); 
theta0 = angle(strain_x + sqrt(-1)*strain_y );


ele = 0;

for k = 1:nelz
    for j = 1:nely
        for i = 1:nelx
            ele = ele + 1;
            epsilon(1,ele) = -strain_xy*abs( cos(theta0+theta(j,i,k) ) );
            epsilon(2,ele) = -strain_xy*abs( sin(theta0+theta(j,i,k) ) );
            epsilon(3,ele) = strain_z;            
            
            if cos( theta0 + theta(j,i,k) ) >= 0
                depsilon(1,ele) = strain_xy*sin( theta0 + theta(j,i,k) );
            else
                depsilon(1,ele) = -strain_xy*sin( theta0 + theta(j,i,k) );
            end
            
            if sin( theta0 + theta(j,i,k) ) >= 0
                depsilon(2,ele) = -strain_xy*cos( theta0 + theta(j,i,k) );
            else
                depsilon(2,ele) = strain_xy*cos( theta0 + theta(j,i,k) );
            end
            
            
        end
    end
end

