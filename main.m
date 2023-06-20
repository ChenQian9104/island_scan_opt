%% =============Island scanning orientation optimization ================%%
%% == objective function: Top layer displacement along Z direction ======%%
%% ========== Finite element solver: Sequential layer-by-layer ==========%%
clc;clear;
nelx = 20; nely = 10; nelz = 5;
nele = nelx*nely*nelz;
nnode = (nelx + 1)*(nely+1)*(nelz+1);
ndofs = 3*nnode; 

dx = 5/1000; dy = 5/1000;  dz = 1.5/1000;

Be = Be_matrix(dx,dy,dz);
% USER-DEFINED MATERIAL PROPERTIES
E = 104E9;           % Young's modulus of solid material
Emin = 1e-9;      % Young's modulus of void-like material
nu = 0.34;         % Poisson's ratio
C = zeros(6,6);
E0 = E/(1+nu)/(1-2*nu);
G = E/(1+nu)/2;
%% For activated elements
C(1,1) = E0*(1-nu);  C(1,2) = E0*nu;    C(1,3) = E0*nu;
C(2,1) = E0*nu;      C(2,2) = E0*(1-nu);C(2,3) = E0*nu;
C(3,1) = E0*nu;      C(3,2) = E0*nu;    C(3,3) = E0*(1-nu);
C(4,4) = G; C(5,5) = G; C(6,6) = G;

edofMat = zeros(nele,8);  % connectivity matrix: [node1 node2 node3 node4 .... node8] - element 1
nodeCor = zeros(nnode, 3); % coordinate of element nodes: [ node1_x node1_y node1_z .... node8_z]
nodeIncrement = (nelx + 1)*(nely + 1);
eleIncrement = nelx*nely;

nstep = nelz + 1;   % simulation steps: nelz layer + cutting off 
for i = 1:nelz 
    element_activation{i} = [ 1: i*eleIncrement ];
end
element_cutting = setdiff( 1:nelx*nely, union( 1:nelx:nelx*nely, 2:nelx:nelx*nely) );
element_activation{nelz + 1} = setdiff( 1:nele, element_cutting);


% define the connectivity matrix for deposition mesh
m=0;
for j = 1:nely 
    for i = 1:nelx 
        m = m+1;       
        node1 = i + (j-1)*(nelx+1);
        node2 = node1+1;
        node3 = node2 + (nelx+1);
        node4 = node1 + (nelx+1);
        node5 = node1 + nodeIncrement;
        node6 = node2 + nodeIncrement;
        node7 = node3 + nodeIncrement;
        node8 = node4 + nodeIncrement;
        edofMat(m,:) = [ node1, node2, node3, node4, node5, node6,node7,node8];
    end
end

for k = 2:nelz
    for i = 1:eleIncrement
        edofMat( i + (k-1)*eleIncrement,:) = edofMat(i + (k-2)*eleIncrement,:) + nodeIncrement; 
    end
end

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

edofMat1 = zeros(nele,24);
for i = 1:8
    edofMat1(:,3*i-2) = 3*edofMat(:,i)-2;
    edofMat1(:,3*i-1) = 3*edofMat(:,i)-1;
    edofMat1(:,3*i) = 3*edofMat(:,i);
end

KE = lk_H8(dx,dy,dz,C);
dlmwrite('stiffness.txt',KE);clear KE;
KE = dlmread('stiffness.txt');

iK = reshape(kron(edofMat1,ones(24,1))',24*24*nele ,1);
jK = reshape(kron(edofMat1,ones(1,24))',24*24*nele ,1);

for i = 1:nstep
    xPhys = Emin*ones(nele,1);
    xPhys( element_activation{i}, 1) = 1.0;
    sK = reshape(KE(:)*xPhys(:)',24*24*nele,1);
    K_layer = sparse(iK,jK,sK); 
    K{i} = (K_layer + K_layer')/2;
end

%% H_act matrix: corresponding to activated elements
H_act{1} = [];
for i = 2:nelz
    H = zeros( ndofs, ndofs);
    activenid = [i*nodeIncrement + 1: (i+1)*nodeIncrement]';
    dofs = [3*activenid(:);3*activenid(:)-1; 3*activenid(:)-2];
    H(dofs, dofs) = 1.0;
    H_act{i} = H;
end 



for i = 1:nelz
    activenid = [nodeIncrement + 1: (i+1)*nodeIncrement]'; 
    activedof{i} = [3*activenid(:); 3*activenid(:)-1; 3*activenid(:)-2]; % DOFs
end
activedof{nstep} = activedof{nelz};


m=0;

F_th = sparse(ndofs,1); % DEFINE THERMAL STRESS
theta = pi/2*ones(nele,1); % scanning angle of each element;
thetaPhys = theta;

%=========================== BOUNDARY CONDITION ==========================%
fixednid = [ 1 : nodeIncrement]';

fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
F = sparse(ndofs,1);
theta = pi/2*ones(nele,1); % scanning angle of each element;


object = [];
sens_1 = []; sens_2 = [];
for iterNum = 1:100
    
    
    %% Linear Elastic Sequential Finite Element Analysis
    [epsilon,depsilon] = Inherent_strain( theta );   % Inherent strain vector
    F_th = sparse(ndofs,1);    
    U = zeros(ndofs,1);
    for i = 1:nstep
        freedofs = activedof{i};

        %% Calculate thermal loading 
        if i == 1
            for ele = 1:nelx*nely
                f_th = thermal_load(Be,C,epsilon(:,i),dx,dy,dz );
                F_th(edofMat1(i,:)',1) = F_th(edofMat1(i,:)',1) + f_th;
 
            end
            U(freedofs,:) = K{i}(freedofs,freedofs)\F_th(freedofs,:);
        else
            continue
        end
        
    end
    
    
    
    
      

    
    compliance = ( F + F_th)'*U;
    object = [object compliance];   
    
    figure(1)
    plot(1:iterNum,object,'b*-','linewidth',2,'markersize',1); pause(0.01)  
    
    dv = sensitivity(Be, C, depsilon, U, edofMat1, nelx, nely,nelz, dx, dy, dz);
    sens_1 = [sens_1, dv(1)];
    snes_2 = [sens_2, dv( nelx*(nely-1) + 1)];
    figure(2)
    plot(1:iterNum,sens_1,'b*-','linewidth',2,'markersize',1);hold on; 
    plot(1:iterNum,sens_2,'r*-','linewidth',2,'markersize',1);hold off; 
    pause(0.01)
    
    if iterNum < 10

        theta = theta - 1*dv;
    elseif iterNum < 20
        theta = theta - 0.5*dv;     
    elseif iterNum < 30
        theta = theta - 0.2*dv;
    else
        theta = theta - 0.02*dv;
    end
    
    theta = mod(theta,pi); 
    
end