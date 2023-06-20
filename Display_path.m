%% Display the scanning path field
u = cos( theta ); v = sin( theta );

figure(4);
u1 = zeros( nely,nelx );
v1 = zeros( nely,nelx  );
[X,Y] = meshgrid(5*0.5:5:5*(nelx-0.5),5*0.5:5:5*(nely-0.5) );

ele = 6*nelx*nely;
for j = 1:10
    for i = 1:20
        ele = ele + 1;
        u1(j,i) = u(ele);
        v1(j,i) = v(ele);
    end
end
q = quiver(X,Y,u1,v1,'MaxHeadSize',10,'AutoScaleFactor',0.4,'AutoScale','on','LineWidth',1.0);
q.Color = 'red';
axis equal;
set(gca,'XLim',[0 100],'XTick',[0:5:100],'YLim',[0 50],'YTick',[0:5:50]);
set(gca,'FontName', 'Times New Roman');
set(gca, 'fontsize', 12);
%set(gca,'fontsize',14);
grid on;



figure(5);
u1 = zeros( nely,nelx );
v1 = zeros( nely,nelx  );
[X,Y] = meshgrid(5*0.5:5:5*(nelx-0.5),5*0.5:5:5*(nely-0.5) );

ele = 8*nelx*nely;
for j = 1:10
    for i = 1:20
        ele = ele + 1;
        u1(j,i) = u(ele);
        v1(j,i) = v(ele);
    end
end
q = quiver(X,Y,u1,v1,'MaxHeadSize',5,'AutoScaleFactor',0.4,'AutoScale','on','LineWidth',1.0);
q.Color = 'red';
axis equal;
set(gca,'XLim',[0 100],'XTick',[0:5:100],'YLim',[0 50],'YTick',[0:5:50]);
set(gca,'FontName', 'Times New Roman');
set(gca, 'fontsize', 12);
grid on;


figure(6);
u1 = zeros( nely,nelx );
v1 = zeros( nely,nelx  );
[X,Y] = meshgrid(5*0.5:5:5*(nelx-0.5),5*0.5:5:5*(nely-0.5) );

ele = 9*nelx*nely;
for j = 1:10
    for i = 1:20
        ele = ele + 1;
        u1(j,i) = u(ele);
        v1(j,i) = v(ele);
    end
end
q = quiver(X,Y,u1,v1,'MaxHeadSize',5,'AutoScaleFactor',0.4,'AutoScale','on','LineWidth',1.0);
q.Color = 'red';
axis equal;
set(gca,'XLim',[0 100],'XTick',[0:5:100],'YLim',[0 50],'YTick',[0:5:50]);
set(gca,'FontName', 'Times New Roman');
set(gca, 'fontsize', 12);
grid on;