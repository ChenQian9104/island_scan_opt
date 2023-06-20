% === GENERATE ELEMENT STIFFNESS MATRIX ===
function [KE] = lk_H8(L1,L2,L3,C)

syms x y z
N1 = (1-2*x/L1)*(1-2*y/L2)*(1-2*z/L3)/8;
N2 = (1+2*x/L1)*(1-2*y/L2)*(1-2*z/L3)/8;
N3 = (1+2*x/L1)*(1+2*y/L2)*(1-2*z/L3)/8;
N4 = (1-2*x/L1)*(1+2*y/L2)*(1-2*z/L3)/8;

N5 = (1-2*x/L1)*(1-2*y/L2)*(1+2*z/L3)/8;
N6 = (1+2*x/L1)*(1-2*y/L2)*(1+2*z/L3)/8;
N7 = (1+2*x/L1)*(1+2*y/L2)*(1+2*z/L3)/8;
N8 = (1-2*x/L1)*(1+2*y/L2)*(1+2*z/L3)/8;

B = zeros(6,24);
B = [diff(N1,x),0,0,diff(N2,x),0,0,diff(N3,x),0,0,diff(N4,x),0,0,...
    diff(N5,x),0,0,diff(N6,x),0,0,diff(N7,x),0,0,diff(N8,x),0,0;...
    0,diff(N1,y),0,0,diff(N2,y),0,0,diff(N3,y),0,0,diff(N4,y),0,...
    0,diff(N5,y),0,0,diff(N6,y),0,0,diff(N7,y),0,0,diff(N8,y),0;...
    0,0,diff(N1,z),0,0,diff(N2,z),0,0,diff(N3,z),0,0,diff(N4,z),...
    0,0,diff(N5,z),0,0,diff(N6,z),0,0,diff(N7,z),0,0,diff(N8,z);...
    diff(N1,y),diff(N1,x),0,diff(N2,y),diff(N2,x),0,diff(N3,y),diff(N3,x),0,diff(N4,y),diff(N4,x),0,...
    diff(N5,y),diff(N5,x),0,diff(N6,y),diff(N6,x),0,diff(N7,y),diff(N7,x),0,diff(N8,y),diff(N8,x),0;
    0,diff(N1,z),diff(N1,y),0,diff(N2,z),diff(N2,y),0,diff(N3,z),diff(N3,y),0,diff(N4,z),diff(N4,y),...
    0,diff(N5,z),diff(N5,y),0,diff(N6,z),diff(N6,y),0,diff(N7,z),diff(N7,y),0,diff(N8,z),diff(N8,y);
    diff(N1,z),0,diff(N1,x),diff(N2,z),0,diff(N2,x),diff(N3,z),0,diff(N3,x),diff(N4,z),0,diff(N4,x),...
    diff(N5,z),0,diff(N5,x),diff(N6,z),0,diff(N6,x),diff(N7,z),0,diff(N7,x),diff(N8,z),0,diff(N8,x);];
    
    
    %K = int( int( N'*N + B'*B,x,-L/2,L/2), y, -L/2,L/2) ;

KE = zeros(24,24);

KE = int( int( int(B'*C*B,x,-L1/2,L1/2), y, -L2/2, L2/2), z, -L3/2, L3/2);




end