function [f_th] = thermal_load(Be,C,Inherent_strain,dx,dy,dz)
f_th = Be'*C*Inherent_strain*dx*dy*dz;
end