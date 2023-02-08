% creation of scatteredinterpolants for Santiago Basin

%% Iquique and SIGAS

clearvars
clc
% load CHI_IQUIQUE.mat
 load CHI_SIGAS.mat
N = numel(spdata);

for i=1:N
    x  = spdata(i).lon;
    y  = spdata(i).lat;
    xy = [x,y];
    xy = unique(xy,'rows','stable');
    
    x = xy(:,1);
    y = xy(:,2);
    
    t  = delaunay(x,y);
    cg  = [mean(x(t),2),mean(y(t),2)];
    OUT = ~inpolygon(cg(:,1),cg(:,2),x,y);
    t(OUT,:)=[];

    n = numel(x);
    spdata(i).lon   = x;
    spdata(i).lat   = y;
    spdata(i).value = spdata(i).value(1)*ones(n,1);
    spdata(i).edge  = (1:n)';
    spdata(i).faces = t;
    spdata(i).F     = scatteredInterpolant(x,y,spdata(i).value,'linear','none');
    i
end
% save CHI_IQUIQUE spdata
 save CHI_SIGAS.mat spdata
%% Cuenca de Santiago
clearvar
Data=xlsread('VS30CuencaEdge.xlsx','DecisionTree');
% Data=xlsread('VS30CuencaEdge.xlsx','LinearRegression');
%  Data=xlsread('VS30CuencaEdge.xlsx','RandomForest');

Edge=xlsread('VS30CuencaEdge.xlsx','BasinEdge');
       
spdata=struct('lon',Data(:,1),'lat',Data(:,2),'value',Data(:,3),'edge',Edge,'faces',[],'F',[]);
t            = delaunay(spdata.lon,spdata.lat);
cg           = [mean(spdata.lon(t),2),mean(spdata.lat(t),2)];
IN           = inpolygon(cg(:,1),cg(:,2),spdata.lon(Edge),spdata.lat(Edge));
spdata.faces = t(IN,:);
spdata.F     = scatteredInterpolant(spdata.lon,spdata.lat,spdata.value,'linear','none');


save SantiagoBasin_DT spdata
% save SantiagoBasin_LR spdata
%  save SantiagoBasin_RF spdata
disp('done')

%% Training Polygons
clearvars
clc
close all
subplot(1,2,1)
N = 2;
spdata(1:N,1)=struct('lon',[],'lat',[],'value',[],'edge',[],'faces',[],'F',[]);

spdata(1).lon=[-1 1 1 -1]';
spdata(1).lat=[0 0 1 1]'+1.1;
spdata(1).value=1;

theta = linspace(0,360,91)';theta(end)=[];
X = cosd(theta);
Y = sind(theta);
Z = theta*0;
[p,t]=mesh_tria([X,Y,Z],0.1,0,[]);
p=p(:,1:2);
% patch('faces',t,'vertices',p,'facecolor','w')
axis equal
edge = zeros(size(theta));
for i=1:numel(edge)
    dist = sum((p-[X(i) Y(i)]).^2,2);
    [val,ind]=min(dist);
    edge(i)=ind;
end
spdata(2).lon=p(:,1);
spdata(2).lat=p(:,2);
spdata(2).value=p(:,1).^2;
spdata(2).edge=edge;
spdata(2).faces=delaunay(p);
spdata(2).F =scatteredInterpolant(spdata(2).lon,spdata(2).lat,spdata(2).value,'linear','none');
patch('faces',spdata(2).faces,'vertices',[spdata(2).lon,spdata(2).lat],'facevertexcdata',spdata(2).value,'facecolor','interp','edgecolor','k')

save Training_spdata spdata

hold on
patch('faces',[1 2 3 4],'vertices',[spdata(1).lon,spdata(1).lat],'facevertexcdata',spdata(1).value*ones(4,1),'facecolor','interp','edgecolor','k')

[xi,yi]=meshgrid(linspace(-1.5,1.5,100),linspace(-1,2.2,100));
IN = inpolygon(xi(:),yi(:),spdata(1).lon,spdata(1).lat);
ui = nan(size(xi));
ui(IN)=spdata(1).value;
vi = spdata(2).F(xi,yi);
title('layer data')


subplot(1,2,2)
hold on
S1=surf(xi,yi,ui);S1.EdgeColor='none';
S2=surf(xi,yi,vi);S2.EdgeColor='none';
axis equal
title('interpolated data')
