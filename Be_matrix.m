function Be = Be_matrix(dx, dy, dz)

a = dx/2; b = dy/2; c = dz/2;
Be = [];
DN=zeros(3,8);

DN(1,1)=-1/a/8;DN(2,1)=-1/b/8;DN(3,1)=-1/c/8;
DN(1,2)= 1/a/8;DN(2,2)=-1/b/8;DN(3,2)=-1/c/8;
DN(1,3)= 1/a/8;DN(2,3)= 1/b/8;DN(3,3)=-1/c/8;
DN(1,4)=-1/a/8;DN(2,4)= 1/b/8;DN(3,4)=-1/c/8;
DN(1,5)=-1/a/8;DN(2,5)=-1/b/8;DN(3,5)= 1/c/8;
DN(1,6)= 1/a/8;DN(2,6)=-1/b/8;DN(3,6)= 1/c/8;
DN(1,7)= 1/a/8;DN(2,7)= 1/b/8;DN(3,7)= 1/c/8;
DN(1,8)=-1/a/8;DN(2,8)= 1/b/8;DN(3,8)= 1/c/8;

for i = 1:8
    B = [ DN(1,i), 0, 0;...
        0, DN(2,i), 0;...
        0, 0, DN(3,i);...
        DN(2,i), DN(1,i), 0;...
        0, DN(3,i), DN(2,i);...
        DN(3,i), 0, DN(1,i);];
    Be = [Be,B];
end
    
    
end


