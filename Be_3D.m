function Be=Be_3D(L)
%% Number of node for each element: nnode
%   2D four node quadrilaterial element: nnode=4;
%   3D eight node hexadral element: nnode=8;
%% Degree of freedom: disDOF
%    2D : disDOF=2*nnode;
%    3D : disDOF=3*nnode;

ele_node=[1 2 3 4 5 6 7 8]';
node_xyz=[0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];

disDOF=3;
strDOF=6;
nnode=8;

intx = 2;   % number of integration points in X direction
inty = 2;   % number of integration points in Y direction
intz = 2;

ie=1;
Be=zeros(strDOF,nnode*disDOF);
for i = 1:intx
    for j = 1:inty
         for k=1:intz
            % Get information about integraton points
            [x,y,z,~] = Gauss_integration(intx,inty,intz,i,j,k);
            
            % Compute shape function N and its local derivative DN at integration points
            [~,DN] = Shape_Quadrangle_8(nnode,x,y,z);
        
            % Compute Jacobi matrix J at integration points
            [~,gDN] = Element_J(ie,nnode,ele_node,node_xyz,DN);
            
            % Compute strain matrix B at integration points
            [B] = Element_B(nnode,disDOF,strDOF,gDN);
            
            Be=Be+B;
    
         end
     end
 end
 Be=Be/8/L;

return;

function [x,y,z,wxyz] = Gauss_integration(intx,inty,intz,i,j,k)

% Gauss  integration constants for 1,2 and 3 points
gxy = zeros(3,3);          % coordinates
wg  = zeros(3,3);          % coefficients of weight

gxy(1,1) =  0.0;
wg(1,1)  =  2.0;
gxy(1,2) = -0.577350269189626;
gxy(2,2) =  0.577350269189626;
wg(1,2)  =  1.0;
wg(2,2)  =  1.0;
gxy(1,3) = -0.774596669241483;
gxy(2,3) =  0.0;
gxy(3,3) =  0.774596669241483;
wg(1,3)  =  0.555555555555556;
wg(2,3)  =  0.888888888888889;
wg(3,3)  =  0.555555555555556;

% Get parameters of integration points
x = gxy(i,intx);
y = gxy(j,inty);
z=  gxy(k,intz);
wxyz = wg(i,intx)*wg(j,inty)*wg(k,intz);

return;


function [N,DN] = Shape_Quadrangle_8(nnode,x,y,z)

N=zeros(1,nnode);
DN=zeros(3,nnode);

N(1)=(1-x)*(1-y)*(1-z)/8;
N(2)=(1+x)*(1-y)*(1-z)/8;
N(3)=(1+x)*(1+y)*(1-z)/8;
N(4)=(1-x)*(1+y)*(1-z)/8;
N(5)=(1-x)*(1-y)*(1+z)/8;
N(6)=(1+x)*(1-y)*(1+z)/8;
N(7)=(1+x)*(1+y)*(1+z)/8;
N(8)=(1-x)*(1+y)*(1+z)/8;

DN(1,1)=-(1-y)*(1-z)/8;DN(2,1)=-(1-x)*(1-z)/8;DN(3,1)=-(1-x)*(1-y)/8;
DN(1,2)=(1-y)*(1-z)/8;DN(2,2)=-(1+x)*(1-z)/8;DN(3,2)=-(1+x)*(1-y)/8;
DN(1,3)=(1+y)*(1-z)/8;DN(2,3)=(1+x)*(1-z)/8;DN(3,3)=-(1+x)*(1+y)/8;
DN(1,4)=-(1+y)*(1-z)/8;DN(2,4)=(1-x)*(1-z)/8;DN(3,4)=-(1-x)*(1+y)/8;
DN(1,5)=-(1-y)*(1+z)/8;DN(2,5)=-(1-x)*(1+z)/8;DN(3,5)=(1-x)*(1-y)/8;
DN(1,6)=(1-y)*(1+z)/8;DN(2,6)=-(1+x)*(1+z)/8;DN(3,6)=(1+x)*(1-y)/8;
DN(1,7)=(1+y)*(1+z)/8;DN(2,7)=(1+x)*(1+z)/8;DN(3,7)=(1+x)*(1+y)/8;
DN(1,8)=-(1+y)*(1+z)/8;DN(2,8)=(1-x)*(1+z)/8;DN(3,8)=(1-x)*(1+y)/8;


return;

%----------------------------------------------
% Compute Jacobi matrix J at integration points
%----------------------------------------------
function [dJ,gDN] = Element_J(ie,nnode,ele_node_dim2,node_xyz_dim2,DN)
    
% Initializaton
J = zeros(3,3);
temp_xyz_dim2 = zeros(nnode,3);

% Matrix of coornidates of nodes in the element
for i = 1:nnode
    temp_xyz_dim2(i,:) = node_xyz_dim2(ele_node_dim2(i,ie),:);
end

% Relevant calculation
J = DN*temp_xyz_dim2;   % the Jacobi matrix
dJ = det(J);            % determinant of the Jacobi matrix            
gDN = J\DN;             % global derivative of N
   
return;

%----------------------------------------------
% Compute strain matrix B at integration points
%----------------------------------------------
function [B] = Element_B(nnode,disDOF,strDOF,gDN)

% Initialization
B = zeros(strDOF,disDOF*nnode);

 for i = 1:nnode
        B(1,(i-1)*3+1) = gDN(1,i);
        B(2,(i-1)*3+2) = gDN(2,i);
        B(3,(i-1)*3+3) = gDN(3,i);
        
        B(4,(i-1)*3+1) = gDN(2,i);
        B(4,(i-1)*3+2) = gDN(1,i);
        
        B(5,(i-1)*3+2) = gDN(3,i);
        B(5,(i-1)*3+3) = gDN(2,i);
        
        B(6,(i-1)*3+1) = gDN(3,i);
        B(6,(i-1)*3+3) = gDN(1,i);

 end
