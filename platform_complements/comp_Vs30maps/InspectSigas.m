clearvars


load CHI_SIGAS.mat
M=zeros(146,1);
hold on
for i=1:146
    M(i)=spdata(i).value(1);
    patch('vertices',[spdata(i).lon,spdata(i).lat],'faces',spdata(i).faces,'facecolor',rand(1,3),'edgecol','none')
    text(mean(spdata(i).lon),mean(spdata(i).lat),num2str(i))
%     spdata(i).lon(end+1)=NaN;
%     spdata(i).lat(end+1)=NaN;
%     spdata(i).value(end+1)=NaN;
    
end
axis equal

uM=unique(M);

spd2(1:numel(uM),1)=struct('lon',[],'lat',[],'value',[],'edge',[],'faces',[],'F',[]);
for i=1:numel(uM)
    IND=M==uM(i);
    spd2(i).lon=vertcat(spdata(IND).lon);
    spd2(i).lat=vertcat(spdata(IND).lat);
end


%%
close all
hold on
patch(spd2(1).lon,spd2(1).lat,rand(1,3))
patch(spd2(2).lon,spd2(2).lat,rand(1,3))
patch(spd2(3).lon,spd2(3).lat,rand(1,3))
patch(spd2(4).lon,spd2(4).lat,rand(1,3))
patch(spd2(5).lon,spd2(5).lat,rand(1,3))
axis equal