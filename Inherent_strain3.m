function [epsilon,depsilon] = Inherent_strain2( theta, island_Mat, nelx, nely,nelz )     % Output the inherent strains of each element
% This transformation is based on the based on the transformation rule on
% the following link:
% https://www.ecourses.ou.edu/cgi-bin/eBook.cgi?doc=&topic=me&chap_sec=08.1&page=theory
epsilon = zeros(6,nelx*nely*nelz); depsilon = zeros(6,nelx*nely*nelz);
[num_island, num_ele] = size( island_Mat );
%strain_x = -0.0145; strain_y = -0.0065; strain_z = 0.012;
strain_x = -0.02; strain_y = -0.01; strain_z = 0.015;
strain_xy = sqrt( strain_x^2 + strain_y^2 ); 
theta0 = angle(strain_x + sqrt(-1)*strain_y );

for m = 1:num_island
    for i = 1:num_ele
        ele = island_Mat(m, i);
        
        epsilon(1,ele) = ( strain_x + strain_y)/2 + (strain_x - strain_y)/2*cos(2*theta(m));
        epsilon(2,ele) = ( strain_x + strain_y)/2 - (strain_x - strain_y)/2*cos(2*theta(m));
        epsilon(3,ele) = strain_z;   
        epsilon(4,ele) = -(strain_x - strain_y)*sin(2*theta(m));
        
        depsilon(1,ele) = -(strain_x - strain_y)*sin( 2*theta(m) );
        depsilon(2,ele) = (strain_x - strain_y)*sin( 2*theta(m) );
        depsilon(4,ele) = -( strain_x - strain_y)*2*cos( 2*theta(m) );
       
        
    end
end


