clearvars

load CHI_USGS.mat
close all
hold on
axis equal
spdata(1:numel(raster))=struct('lon',[],'lat',[],'value',[],'edge',[],'faces',[],'F',[]);
idx=[];
raster([1,19,13,21])=[];
for i=1:numel(raster)
    r     = raster(i);
    img   = double(imread(r.tif));
    lon   = linspace(r.box(1,1),r.box(2,1),r.nlon);
    lat   = linspace(r.box(2,2),r.box(1,2),r.nlat);
    [LON,LAT] = meshgrid(lon,lat);
    LON = LON(:);
    LAT = LAT(:);
    VAL = img(:);
    if numel(lon)>=3 && numel(lat)>=3
        t  = delaunay(LON,LAT);
        patch('Faces',t,...
            'Vertices',[LON,LAT],...
            'facevertexcdata',VAL,...
            'facecol','interp',...
            'edgecol','none',...
            'facealpha',0.6,...
            'tag','layer')
        text(mean(LON),mean(LAT),num2str(i))
        
        % populates spdata
        spdata(i).lon=LON;
        spdata(i).lat=LAT;
        spdata(i).value=VAL;
        
        lonE = lon([1,end,end,1])';
        latE = lat([1,1,end,end])';
        [~,b]=intersect([LON,LAT],[lonE,latE],'rows','stable');
        b = b([1 2 4 3]);
        spdata(i).edge=b;
        spdata(i).faces=t;
        spdata(i).F = scatteredInterpolant(LON,LAT,VAL,'linear','none');
        
        edge = b;
        plot(spdata(i).lon(edge),spdata(i).lat(edge),'r-')
        plot(spdata(i).lon(edge(1)),spdata(i).lat(edge(1)),'ro')
    else
        idx=[idx;i];
        
    end
    i
end
idx;
spdata(idx)=[];
save CHI_USGS_develop spdata
%%
load CHI_USGS.mat
load CHI_USGS_develop.mat
raster([1,19,13,21])=[];
BOX = vertcat(raster.box);
box  = [min(BOX(:,1)),min(BOX(:,2));
        max(BOX(:,1)),max(BOX(:,2))];
nx = round(diff(box(:,1))/0.0084);
ny = round(diff(box(:,2))/0.0084);
x  = linspace(box(1,1),box(2,1),nx);
y  = linspace(box(1,2),box(2,2),ny);
[X,Y]=meshgrid(x,y);
X = X(:);
Y = Y(:);
V = nan(length(X),length(spdata));

for i=1:length(spdata)
    if ~isempty(spdata(i).F)
        edge = spdata(i).edge;
        px = spdata(i).lon(edge);
        py = spdata(i).lat(edge);
        IN = inpolygon(X,Y,px,py);
        V(IN,i)=spdata(i).F(X(IN),Y(IN));
    end
end

V = mean(V,2,'omitnan');
OUT = isnan(V);

X(OUT)=[];
Y(OUT)=[];
V(OUT)=[];
plot(X,Y,'.');
axis equal

save myTEMPDATA X Y V
%%
clearvars
clc
load myTEMPDATA.mat
close all,hold on,axis equal
plot(X,Y,'.','ButtonDownFcn',{@getsitepointer,X,Y});
edge=[     1453574
     1453373
     1482229
     1482020
     1335724
     1335206
     1068804
     1068358
      998008
      997653
      954873
      954655
      888805
      888706
      668161
      668004
      343075
      342780
      460823
      460727
      403652
      403539
      147883
      147478
        4985
        4395
   879
     1
     2337583
     2338461
     1228483
     1229073
     1047582
     1047987
     1044780
     1044853
     1307008
     1307251
     1260349
     1260474
     1469144
     1469327
     1575163
     1575248
     1655578
     1655740
     1586922
     1586999
     1643039
     1643366
     2052709
     2053123
     2276323
     2276931
     2021159
     2021397
     1892825
     1892997];
 
 
 t  = delaunay(X,Y);
 cg  = [mean(X(t),2),mean(Y(t),2)];
 OUT = ~inpolygon(cg(:,1),cg(:,2),X(edge),Y(edge));
 t(OUT,:)=[];
 
 
 spdata.lon   = X;
 spdata.lat   = Y;
 spdata.value = V;
 spdata.edge  = edge;
 spdata.faces = t;
 spdata.F     = scatteredInterpolant(X,Y,V,'linear','none');



