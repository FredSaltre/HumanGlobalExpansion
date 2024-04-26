program DataSource_FigS5_TableS2

% we are interested in column 4,5 and 6
%%=========================================================================
%% CALCULATING THE LENGTH & SPEED OF ALL TRAJECTORIES + COMPARING THE BEST ONE WITH OTHERS
%%=========================================================================
%% AFRICA-BERINGIA
%% the two best trajectoris are here #3 and #5
clear all, close all, clc
pth1=load('Trajectory_Africa-Beringia_(path#1).txt');pth2=load('Trajectory_Africa-Beringia_(path#2).txt');pth3=load('Trajectory_Africa-Beringia_(path#3).txt');
pth4=load('Trajectory_Africa-Beringia_(path#4).txt');pth5=load('Trajectory_Africa-Beringia_(path#5).txt');pth6=load('Trajectory_Africa-Beringia_(path#6).txt');

% PATH 1
path1.dat=sortrows(pth1(:,4:6),3,'descend');
path1.dif=abs(diff(path1.dat(:,3)));% time difference beteeen grid cell
path1.cum=sum(path1.dat(:,3));%total time along the trajectory
path1.dist=nan(length(path1.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path1.dat)-1);path1.dist(i) = deg2km(distance(path1.dat(i,1),path1.dat(i,2),path1.dat(i+1,1),path1.dat(i+1,2)));end;
path1.speed=path1.dist./path1.dif;

outpath1=[ones(length(pth1(:,1)),1),pth1(:,4:6),[0;path1.dist],[0;path1.speed]];


% PATH 2
path2.dat=sortrows(pth2(:,4:6),3,'descend');
path2.dif=abs(diff(path2.dat(:,3)));% time difference beteeen grid cell
path2.cum=sum(path2.dat(:,3));%total time along the trajectory
path2.dist=nan(length(path2.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path2.dat)-1);path2.dist(i) = deg2km(distance(path2.dat(i,1),path2.dat(i,2),path2.dat(i+1,1),path2.dat(i+1,2)));end;
path2.speed=path2.dist./path2.dif;%figure,histogram(path1.speed,100)

outpath2=[repmat(2,length(pth2(:,1)),1),pth2(:,4:6),[0;path2.dist],[0;path2.speed]];


% PATH 3
path3.dat=sortrows(pth3(:,4:6),3,'descend');
path3.dif=abs(diff(path3.dat(:,3)));% time difference beteeen grid cell
path3.cum=sum(path3.dat(:,3));%total time along the trajectory
path3.dist=nan(length(path3.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path3.dat)-1);path3.dist(i) = deg2km(distance(path3.dat(i,1),path3.dat(i,2),path3.dat(i+1,1),path3.dat(i+1,2)));end;
path3.speed=path3.dist./path3.dif;%figure,histogram(path1.speed,100)

outpath3=[repmat(3,length(pth3(:,1)),1),pth3(:,4:6),[0;path3.dist],[0;path3.speed]];

% PATH 4
path4.dat=sortrows(pth4(:,4:6),3,'descend');
path4.dif=abs(diff(path4.dat(:,3)));% time difference beteeen grid cell
path4.cum=sum(path4.dat(:,3));%total time along the trajectory
path4.dist=nan(length(path4.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path4.dat)-1);path4.dist(i) = deg2km(distance(path4.dat(i,1),path4.dat(i,2),path4.dat(i+1,1),path4.dat(i+1,2)));end;
path4.speed=path4.dist./path4.dif;%figure,histogram(path4.speed,100)

outpath4=[repmat(4,length(pth4(:,1)),1),pth4(:,4:6),[0;path4.dist],[0;path4.speed]];

% PATH 5
path5.dat=sortrows(pth5(:,4:6),3,'descend');
path5.dif=abs(diff(path5.dat(:,3)));% time difference beteeen grid cell
path5.cum=sum(path5.dat(:,3));%total time along the trajectory
path5.dist=nan(length(path5.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path5.dat)-1);path5.dist(i) = deg2km(distance(path5.dat(i,1),path5.dat(i,2),path5.dat(i+1,1),path5.dat(i+1,2)));end;
path5.speed=path5.dist./path5.dif;%figure,histogram(path5.speed,100)

outpath5=[repmat(5,length(pth5(:,1)),1),pth5(:,4:6),[0;path5.dist],[0;path5.speed]];

% PATH 6
path6.dat=sortrows(pth6(:,4:6),3,'descend');
path6.dif=abs(diff(path6.dat(:,3)));% time difference beteeen grid cell
path6.cum=sum(path6.dat(:,3));%total time along the trajectory
path6.dist=nan(length(path6.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path6.dat)-1);path6.dist(i) = deg2km(distance(path6.dat(i,1),path6.dat(i,2),path6.dat(i+1,1),path6.dat(i+1,2)));end;
path6.speed=path6.dist./path6.dif;%figure,histogram(path6.speed,100)

outpath6=[repmat(6,length(pth6(:,1)),1),pth6(:,4:6),[0;path6.dist],[0;path6.speed]];

out=[outpath1;outpath2;outpath3;outpath4;outpath5;outpath6];nout=length(out(:,1));
gen=load('Humancolonisation_FSTdata_v2.txt');out=[out,nan(nout,1)];%load genetic data
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg')
%coast = load('coastlines');seaclr1=hex2rgb('F0F8FF');%convert hexadecimal color into RGB triplet
%figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
%                    geoshow(coast.coastlat,coast.coastlon,'DisplayType','polygon','Facecolor','w'),
%                    geoshow(flipud(coast.coastlat),flipud(coast.coastlon),'DisplayType','polygon','FaceColor',seaclr1),
%                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
%                    geoshow(land, 'FaceColor', [0.8 0.8 0.8]);hold on,alpha(0.7);
%                    hold on,surfm(yg+0.5,xg2-0.5,gen);
for i=1:nout;
    vlat=abs(yg-out(i,2));%idlat=find(vlat==min(vlat));
    vlon=abs(xg2-out(i,3));%idlon=find(vlon==min(vlon));
    out(i,end)=gen(vlat==min(vlat),vlon==min(vlon));
    clear vlat vlon;% idlat idlon;
end;

% Now create the table
FigSD = table(out(:,1), out(:,2), out(:,3), out(:,4), out(:,5), out(:,6), out(:,7),'VariableNames', {'Trajectory', 'Latitude', 'Longitude', 'Date', 'Length','Speed', 'FST'});
% Save the table to a CSV file
writetable(FigSD, 'SourceData(Africa_Beringia).csv');

%saving the outcome
AfrBer.path1=path1;AfrBer.path2=path2;AfrBer.path3=path3;AfrBer.path4=path4;AfrBer.path5=path5;AfrBer.path6=path6;
save('Africa_Beringia.mat', '-struct', 'AfrBer');


%%=========================================================================
%% MEXICO-CHILE
%% the two best trajectoris are here #1 and #3
clear all, close all, clc,       
pth1=load('Trajectory_Mexico-Chile_(path#1).txt');pth2=load('Trajectory_Mexico-Chile_(path#2).txt');pth3=load('Trajectory_Mexico-Chile_(path#3).txt');
pth4=load('Trajectory_Mexico-Chile_(path#4).txt');pth5=load('Trajectory_Mexico-Chile_(path#5).txt');pth6=load('Trajectory_Mexico-Chile_(path#6).txt');
pth7=load('Trajectory_Mexico-Chile_(path#7).txt');

% PATH 1
path1.dat=sortrows(pth1(:,4:6),3,'descend');
path1.dif=abs(diff(path1.dat(:,3)));% time difference beteeen grid cell
path1.cum=sum(path1.dat(:,3));%total time along the trajectory
path1.dist=nan(length(path1.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path1.dat)-1);path1.dist(i) = deg2km(distance(path1.dat(i,1),path1.dat(i,2),path1.dat(i+1,1),path1.dat(i+1,2)));end;
path1.speed=path1.dist./path1.dif;%figure,histogram(path1.speed,100)

outpath1=[ones(length(pth1(:,1)),1),pth1(:,4:6),[0;path1.dist],[0;path1.speed]];

% PATH 2
path2.dat=sortrows(pth2(:,4:6),3,'descend');
path2.dif=abs(diff(path2.dat(:,3)));% time difference beteeen grid cell
path2.cum=sum(path2.dat(:,3));%total time along the trajectory
path2.dist=nan(length(path2.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path2.dat)-1);path2.dist(i) = deg2km(distance(path2.dat(i,1),path2.dat(i,2),path2.dat(i+1,1),path2.dat(i+1,2)));end;
path2.speed=path2.dist./path2.dif;%figure,histogram(path1.speed,100)

outpath2=[repmat(2,length(pth2(:,1)),1),pth2(:,4:6),[0;path2.dist],[0;path2.speed]];

% PATH 3
path3.dat=sortrows(pth3(:,4:6),3,'descend');
path3.dif=abs(diff(path3.dat(:,3)));% time difference beteeen grid cell
path3.cum=sum(path3.dat(:,3));%total time along the trajectory
path3.dist=nan(length(path3.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path3.dat)-1);path3.dist(i) = deg2km(distance(path3.dat(i,1),path3.dat(i,2),path3.dat(i+1,1),path3.dat(i+1,2)));end;
path3.speed=path3.dist./path3.dif;%figure,histogram(path1.speed,100)

outpath3=[repmat(3,length(pth3(:,1)),1),pth3(:,4:6),[0;path3.dist],[0;path3.speed]];

% PATH 4
path4.dat=sortrows(pth4(:,4:6),3,'descend');
path4.dif=abs(diff(path4.dat(:,3)));% time difference beteeen grid cell
path4.cum=sum(path4.dat(:,3));%total time along the trajectory
path4.dist=nan(length(path4.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path4.dat)-1);path4.dist(i) = deg2km(distance(path4.dat(i,1),path4.dat(i,2),path4.dat(i+1,1),path4.dat(i+1,2)));end;
path4.speed=path4.dist./path4.dif;%figure,histogram(path4.speed,100)

outpath4=[repmat(4,length(pth4(:,1)),1),pth4(:,4:6),[0;path4.dist],[0;path4.speed]];

% PATH 5
path5.dat=sortrows(pth5(:,4:6),3,'descend');
path5.dif=abs(diff(path5.dat(:,3)));% time difference beteeen grid cell
path5.cum=sum(path5.dat(:,3));%total time along the trajectory
path5.dist=nan(length(path5.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path5.dat)-1);path5.dist(i) = deg2km(distance(path5.dat(i,1),path5.dat(i,2),path5.dat(i+1,1),path5.dat(i+1,2)));end;
path5.speed=path5.dist./path5.dif;%figure,histogram(path5.speed,100)

outpath5=[repmat(5,length(pth5(:,1)),1),pth5(:,4:6),[0;path5.dist],[0;path5.speed]];

% PATH 6
path6.dat=sortrows(pth6(:,4:6),3,'descend');
path6.dif=abs(diff(path6.dat(:,3)));% time difference beteeen grid cell
path6.cum=sum(path6.dat(:,3));%total time along the trajectory
path6.dist=nan(length(path6.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path6.dat)-1);path6.dist(i) = deg2km(distance(path6.dat(i,1),path6.dat(i,2),path6.dat(i+1,1),path6.dat(i+1,2)));end;
path6.speed=path6.dist./path6.dif;%figure,histogram(path6.speed,100)

outpath6=[repmat(6,length(pth6(:,1)),1),pth6(:,4:6),[0;path6.dist],[0;path6.speed]];

% PATH 7
path7.dat=sortrows(pth7(:,4:6),3,'descend');
path7.dif=abs(diff(path7.dat(:,3)));% time difference beteeen grid cell
path7.cum=sum(path7.dat(:,3));%total time along the trajectory
path7.dist=nan(length(path7.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path7.dat)-1);path7.dist(i) = deg2km(distance(path7.dat(i,1),path7.dat(i,2),path7.dat(i+1,1),path7.dat(i+1,2)));end;
path7.speed=path7.dist./path7.dif;%figure,histogram(path7.speed,100)

outpath7=[repmat(7,length(pth7(:,1)),1),pth7(:,4:6),[0;path7.dist],[0;path7.speed]];

out=[outpath1;outpath2;outpath3;outpath4;outpath5;outpath6;outpath7];nout=length(out(:,1));
gen=load('Humancolonisation_FSTdata_v2.txt');out=[out,nan(nout,1)];%load genetic data
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg')
%coast = load('coastlines');seaclr1=hex2rgb('F0F8FF');%convert hexadecimal color into RGB triplet
%figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
%                    geoshow(coast.coastlat,coast.coastlon,'DisplayType','polygon','Facecolor','w'),
%                    geoshow(flipud(coast.coastlat),flipud(coast.coastlon),'DisplayType','polygon','FaceColor',seaclr1),
%                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
%                    geoshow(land, 'FaceColor', [0.8 0.8 0.8]);hold on,alpha(0.7);
%                    hold on,surfm(yg+0.5,xg2-0.5,gen);
for i=1:nout;
    vlat=abs(yg-out(i,2));%idlat=find(vlat==min(vlat));
    vlon=abs(xg2-out(i,3));%idlon=find(vlon==min(vlon));
    out(i,end)=gen(vlat==min(vlat),vlon==min(vlon));
    clear vlat vlon;% idlat idlon;
end;

% Now create the table
FigSD = table(out(:,1), out(:,2), out(:,3), out(:,4), out(:,5), out(:,6), out(:,7),'VariableNames', {'Trajectory', 'Latitude', 'Longitude', 'Date', 'Length', 'Speed', 'FST'});
% Save the table to a CSV file
writetable(FigSD, 'SourceData(MexicoChile).csv');


%saving the outcome
MexChil.path1=path1;MexChil.path2=path2;MexChil.path3=path3;MexChil.path4=path4;MexChil.path5=path5;MexChil.path6=path6;MexChil.path7=path7;
save('MexicoChile.mat', '-struct', 'MexChil');

                    
%%=========================================================================
%% AFRICA-FINLAND
%% the two best trajectoris are here #1 and #2
clear all, close all, clc,
pth1=load('Trajectory_Africa-Finland_(path#1).txt');pth2=load('Trajectory_Africa-Finland_(path#2).txt');pth3=load('Trajectory_Africa-Finland_(path#3).txt');
pth4=load('Trajectory_Africa-Finland_(path#4).txt');pth5=load('Trajectory_Africa-Finland_(path#5).txt');pth6=load('Trajectory_Africa-Finland_(path#6).txt');
pth7=load('Trajectory_Africa-Finland_(path#7).txt');

% PATH 1
path1.dat=sortrows(pth1(:,4:6),3,'descend');
path1.dif=abs(diff(path1.dat(:,3)));% time difference beteeen grid cell
path1.cum=sum(path1.dat(:,3));%total time along the trajectory
path1.dist=nan(length(path1.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path1.dat)-1);path1.dist(i) = deg2km(distance(path1.dat(i,1),path1.dat(i,2),path1.dat(i+1,1),path1.dat(i+1,2)));end;
path1.speed=path1.dist./path1.dif;%figure,histogram(path1.speed,100)

outpath1=[ones(length(pth1(:,1)),1),pth1(:,4:6),[0;path1.dist],[0;path1.speed]];

% PATH 2
path2.dat=sortrows(pth2(:,4:6),3,'descend');
path2.dif=abs(diff(path2.dat(:,3)));% time difference beteeen grid cell
path2.cum=sum(path2.dat(:,3));%total time along the trajectory
path2.dist=nan(length(path2.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path2.dat)-1);path2.dist(i) = deg2km(distance(path2.dat(i,1),path2.dat(i,2),path2.dat(i+1,1),path2.dat(i+1,2)));end;
path2.speed=path2.dist./path2.dif;%figure,histogram(path1.speed,100)

outpath2=[repmat(2,length(pth2(:,1)),1),pth2(:,4:6),[0;path2.dist],[0;path2.speed]];

% PATH 3
path3.dat=sortrows(pth3(:,4:6),3,'descend');
path3.dif=abs(diff(path3.dat(:,3)));% time difference beteeen grid cell
path3.cum=sum(path3.dat(:,3));%total time along the trajectory
path3.dist=nan(length(path3.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path3.dat)-1);path3.dist(i) = deg2km(distance(path3.dat(i,1),path3.dat(i,2),path3.dat(i+1,1),path3.dat(i+1,2)));end;
path3.speed=path3.dist./path3.dif;%figure,histogram(path1.speed,100)

outpath3=[repmat(3,length(pth3(:,1)),1),pth3(:,4:6),[0;path3.dist],[0;path3.speed]];

% PATH 4
path4.dat=sortrows(pth4(:,4:6),3,'descend');
path4.dif=abs(diff(path4.dat(:,3)));% time difference beteeen grid cell
path4.cum=sum(path4.dat(:,3));%total time along the trajectory
path4.dist=nan(length(path4.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path4.dat)-1);path4.dist(i) = deg2km(distance(path4.dat(i,1),path4.dat(i,2),path4.dat(i+1,1),path4.dat(i+1,2)));end;
path4.speed=path4.dist./path4.dif;%figure,histogram(path4.speed,100)

outpath4=[repmat(4,length(pth4(:,1)),1),pth4(:,4:6),[0;path4.dist],[0;path4.speed]];

% PATH 5
path5.dat=sortrows(pth5(:,4:6),3,'descend');
path5.dif=abs(diff(path5.dat(:,3)));% time difference beteeen grid cell
path5.cum=sum(path5.dat(:,3));%total time along the trajectory
path5.dist=nan(length(path5.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path5.dat)-1);path5.dist(i) = deg2km(distance(path5.dat(i,1),path5.dat(i,2),path5.dat(i+1,1),path5.dat(i+1,2)));end;
path5.speed=path5.dist./path5.dif;%figure,histogram(path5.speed,100)

outpath5=[repmat(5,length(pth5(:,1)),1),pth5(:,4:6),[0;path5.dist],[0;path5.speed]];

% PATH 6
path6.dat=sortrows(pth6(:,4:6),3,'descend');
path6.dif=abs(diff(path6.dat(:,3)));% time difference beteeen grid cell
path6.cum=sum(path6.dat(:,3));%total time along the trajectory
path6.dist=nan(length(path6.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path6.dat)-1);path6.dist(i) = deg2km(distance(path6.dat(i,1),path6.dat(i,2),path6.dat(i+1,1),path6.dat(i+1,2)));end;
path6.speed=path6.dist./path6.dif;%figure,histogram(path6.speed,100)

outpath6=[repmat(6,length(pth6(:,1)),1),pth6(:,4:6),[0;path6.dist],[0;path6.speed]];

% PATH 7
path7.dat=sortrows(pth7(:,4:6),3,'descend');
path7.dif=abs(diff(path7.dat(:,3)));% time difference beteeen grid cell
path7.cum=sum(path7.dat(:,3));%total time along the trajectory
path7.dist=nan(length(path7.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path7.dat)-1);path7.dist(i) = deg2km(distance(path7.dat(i,1),path7.dat(i,2),path7.dat(i+1,1),path7.dat(i+1,2)));end;
path7.speed=path7.dist./path7.dif;%figure,histogram(path7.speed,100)

outpath7=[repmat(7,length(pth7(:,1)),1),pth7(:,4:6),[0;path7.dist],[0;path7.speed]];

out=[outpath1;outpath2;outpath3;outpath4;outpath5;outpath6;outpath7];nout=length(out(:,1));
gen=load('Humancolonisation_FSTdata_v2.txt');out=[out,nan(nout,1)];%load genetic data
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg')
%coast = load('coastlines');seaclr1=hex2rgb('F0F8FF');%convert hexadecimal color into RGB triplet
%figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
%                    geoshow(coast.coastlat,coast.coastlon,'DisplayType','polygon','Facecolor','w'),
%                    geoshow(flipud(coast.coastlat),flipud(coast.coastlon),'DisplayType','polygon','FaceColor',seaclr1),
%                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
%                    geoshow(land, 'FaceColor', [0.8 0.8 0.8]);hold on,alpha(0.7);
%                    hold on,surfm(yg+0.5,xg2-0.5,gen);
for i=1:nout;
    vlat=abs(yg-out(i,2));%idlat=find(vlat==min(vlat));
    vlon=abs(xg2-out(i,3));%idlon=find(vlon==min(vlon));
    out(i,end)=gen(vlat==min(vlat),vlon==min(vlon));
    clear vlat vlon;% idlat idlon;
end;

% Now create the table
FigSD = table(out(:,1), out(:,2), out(:,3), out(:,4), out(:,5), out(:,6), out(:,7),'VariableNames', {'Trajectory', 'Latitude', 'Longitude', 'Date', 'Length', 'Speed', 'FST'});
% Save the table to a CSV file
writetable(FigSD, 'SourceData(AfricaFinland).csv');

%saving the outcome
AfrFin.path1=path1;AfrFin.path2=path2;AfrFin.path3=path3;AfrFin.path4=path4;AfrFin.path5=path5;AfrFin.path6=path6;AfrFin.path7=path7;
save('AfricaFinland.mat', '-struct', 'AfrFin');


%%=========================================================================
%% AFRICA-JAPAN
%% the two best trajectoris are here #1 and #2
clear all, close all, clc,

pth1=load('Trajectory_Africa-Japan_(path#1).txt');pth2=load('Trajectory_Africa-Japan_(path#2).txt');pth3=load('Trajectory_Africa-Japan_(path#3).txt');
pth4=load('Trajectory_Africa-Japan_(path#4).txt');pth5=load('Trajectory_Africa-Japan_(path#5).txt');pth6=load('Trajectory_Africa-Japan_(path#6).txt');
pth7=load('Trajectory_Africa-Japan_(path#7).txt');pth8=load('Trajectory_Africa-Japan_(path#8).txt');

% PATH 1
path1.dat=sortrows(pth1(:,4:6),3,'descend');
path1.dif=abs(diff(path1.dat(:,3)));% time difference beteeen grid cell
path1.cum=sum(path1.dat(:,3));%total time along the trajectory
path1.dist=nan(length(path1.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path1.dat)-1);path1.dist(i) = deg2km(distance(path1.dat(i,1),path1.dat(i,2),path1.dat(i+1,1),path1.dat(i+1,2)));end;
path1.speed=path1.dist./path1.dif;%figure,histogram(path1.speed,100)

outpath1=[ones(length(pth1(:,1)),1),pth1(:,4:6),[0;path1.dist],[0;path1.speed]];

% PATH 2
path2.dat=sortrows(pth2(:,4:6),3,'descend');
path2.dif=abs(diff(path2.dat(:,3)));% time difference beteeen grid cell
path2.cum=sum(path2.dat(:,3));%total time along the trajectory
path2.dist=nan(length(path2.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path2.dat)-1);path2.dist(i) = deg2km(distance(path2.dat(i,1),path2.dat(i,2),path2.dat(i+1,1),path2.dat(i+1,2)));end;
path2.speed=path2.dist./path2.dif;%figure,histogram(path1.speed,100)

outpath2=[repmat(2,length(pth2(:,1)),1),pth2(:,4:6),[0;path2.dist],[0;path2.speed]];

% PATH 3
path3.dat=sortrows(pth3(:,4:6),3,'descend');
path3.dif=abs(diff(path3.dat(:,3)));% time difference beteeen grid cell
path3.cum=sum(path3.dat(:,3));%total time along the trajectory
path3.dist=nan(length(path3.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path3.dat)-1);path3.dist(i) = deg2km(distance(path3.dat(i,1),path3.dat(i,2),path3.dat(i+1,1),path3.dat(i+1,2)));end;
path3.speed=path3.dist./path3.dif;%figure,histogram(path1.speed,100)

outpath3=[repmat(3,length(pth3(:,1)),1),pth3(:,4:6),[0;path3.dist],[0;path3.speed]];

% PATH 4
path4.dat=sortrows(pth4(:,4:6),3,'descend');
path4.dif=abs(diff(path4.dat(:,3)));% time difference beteeen grid cell
path4.cum=sum(path4.dat(:,3));%total time along the trajectory
path4.dist=nan(length(path4.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path4.dat)-1);path4.dist(i) = deg2km(distance(path4.dat(i,1),path4.dat(i,2),path4.dat(i+1,1),path4.dat(i+1,2)));end;
path4.speed=path4.dist./path4.dif;%figure,histogram(path4.speed,100)

outpath4=[repmat(4,length(pth4(:,1)),1),pth4(:,4:6),[0;path4.dist],[0;path4.speed]];

% PATH 5
path5.dat=sortrows(pth5(:,4:6),3,'descend');
path5.dif=abs(diff(path5.dat(:,3)));% time difference beteeen grid cell
path5.cum=sum(path5.dat(:,3));%total time along the trajectory
path5.dist=nan(length(path5.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path5.dat)-1);path5.dist(i) = deg2km(distance(path5.dat(i,1),path5.dat(i,2),path5.dat(i+1,1),path5.dat(i+1,2)));end;
path5.speed=path5.dist./path5.dif;%figure,histogram(path5.speed,100)

outpath5=[repmat(5,length(pth5(:,1)),1),pth5(:,4:6),[0;path5.dist],[0;path5.speed]];

% PATH 6
path6.dat=sortrows(pth6(:,4:6),3,'descend');
path6.dif=abs(diff(path6.dat(:,3)));% time difference beteeen grid cell
path6.cum=sum(path6.dat(:,3));%total time along the trajectory
path6.dist=nan(length(path6.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path6.dat)-1);path6.dist(i) = deg2km(distance(path6.dat(i,1),path6.dat(i,2),path6.dat(i+1,1),path6.dat(i+1,2)));end;
path6.speed=path6.dist./path6.dif;%figure,histogram(path6.speed,100)

outpath6=[repmat(6,length(pth6(:,1)),1),pth6(:,4:6),[0;path6.dist],[0;path6.speed]];

% PATH 7
path7.dat=sortrows(pth7(:,4:6),3,'descend');
path7.dif=abs(diff(path7.dat(:,3)));% time difference beteeen grid cell
path7.cum=sum(path7.dat(:,3));%total time along the trajectory
path7.dist=nan(length(path7.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path7.dat)-1);path7.dist(i) = deg2km(distance(path7.dat(i,1),path7.dat(i,2),path7.dat(i+1,1),path7.dat(i+1,2)));end;
path7.speed=path7.dist./path7.dif;%figure,histogram(path7.speed,100)

outpath7=[repmat(7,length(pth7(:,1)),1),pth7(:,4:6),[0;path7.dist],[0;path7.speed]];

% PATH 8
path8.dat=sortrows(pth8(:,4:6),3,'descend');
path8.dif=abs(diff(path8.dat(:,3)));% time difference beteeen grid cell
path8.cum=sum(path8.dat(:,3));%total time along the trajectory
path8.dist=nan(length(path8.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path8.dat)-1);path8.dist(i) = deg2km(distance(path8.dat(i,1),path8.dat(i,2),path8.dat(i+1,1),path8.dat(i+1,2)));end;
path8.speed=path8.dist./path8.dif;%figure,histogram(path8.speed,100)

outpath8=[repmat(8,length(pth8(:,1)),1),pth8(:,4:6),[0;path8.dist],[0;path8.speed]];

out=[outpath1;outpath2;outpath3;outpath4;outpath5;outpath6;outpath7;outpath8];nout=length(out(:,1));
gen=load('Humancolonisation_FSTdata_v2.txt');out=[out,nan(nout,1)];%load genetic data
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg')
%coast = load('coastlines');seaclr1=hex2rgb('F0F8FF');%convert hexadecimal color into RGB triplet
%figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
%                    geoshow(coast.coastlat,coast.coastlon,'DisplayType','polygon','Facecolor','w'),
%                    geoshow(flipud(coast.coastlat),flipud(coast.coastlon),'DisplayType','polygon','FaceColor',seaclr1),
%                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
%                    geoshow(land, 'FaceColor', [0.8 0.8 0.8]);hold on,alpha(0.7);
%                    hold on,surfm(yg+0.5,xg2-0.5,gen);
for i=1:nout;
    vlat=abs(yg-out(i,2));%idlat=find(vlat==min(vlat));
    vlon=abs(xg2-out(i,3));%idlon=find(vlon==min(vlon));
    out(i,end)=gen(vlat==min(vlat),vlon==min(vlon));
    clear vlat vlon;% idlat idlon;
end;

% Now create the table
FigSD = table(out(:,1), out(:,2), out(:,3), out(:,4), out(:,5), out(:,6), out(:,7),'VariableNames', {'Trajectory', 'Latitude', 'Longitude', 'Date', 'Length', 'Speed', 'FST'});
% Save the table to a CSV file
writetable(FigSD, 'SourceData(AfricaJapan).csv');

%saving the outcome
AfrJap.path1=path1;AfrJap.path2=path2;AfrJap.path3=path3;AfrJap.path4=path4;AfrJap.path5=path5;AfrJap.path6=path6;AfrJap.path7=path7;AfrJap.path8=path8;
save('AfricaJapan.mat', '-struct', 'AfrJap');

                    
%%=========================================================================
%% AFRICA-PORTUGAL
%% the two best trajectoris are here #1 and #2
clear all, close all, clc,
pth1=load('Trajectory_Africa-Portugal_(path#1).txt');pth2=load('Trajectory_Africa-Portugal_(path#2).txt');pth3=load('Trajectory_Africa-Portugal_(path#3).txt');
pth4=load('Trajectory_Africa-Portugal_(path#4).txt');pth5=load('Trajectory_Africa-Portugal_(path#5).txt');pth6=load('Trajectory_Africa-Portugal_(path#6).txt');

% PATH 1
path1.dat=sortrows(pth1(:,4:6),3,'descend');
path1.dif=abs(diff(path1.dat(:,3)));% time difference beteeen grid cell
path1.cum=sum(path1.dat(:,3));%total time along the trajectory
path1.dist=nan(length(path1.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path1.dat)-1);path1.dist(i) = deg2km(distance(path1.dat(i,1),path1.dat(i,2),path1.dat(i+1,1),path1.dat(i+1,2)));end;
path1.speed=path1.dist./path1.dif;%figure,histogram(path1.speed,100)

outpath1=[ones(length(pth1(:,1)),1),pth1(:,4:6),[0;path1.dist],[0;path1.speed]];

% PATH 2
path2.dat=sortrows(pth2(:,4:6),3,'descend');
path2.dif=abs(diff(path2.dat(:,3)));% time difference beteeen grid cell
path2.cum=sum(path2.dat(:,3));%total time along the trajectory
path2.dist=nan(length(path2.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path2.dat)-1);path2.dist(i) = deg2km(distance(path2.dat(i,1),path2.dat(i,2),path2.dat(i+1,1),path2.dat(i+1,2)));end;
path2.speed=path2.dist./path2.dif;%figure,histogram(path1.speed,100)

outpath2=[repmat(2,length(pth2(:,1)),1),pth2(:,4:6),[0;path2.dist],[0;path2.speed]];

% PATH 3
path3.dat=sortrows(pth3(:,4:6),3,'descend');
path3.dif=abs(diff(path3.dat(:,3)));% time difference beteeen grid cell
path3.cum=sum(path3.dat(:,3));%total time along the trajectory
path3.dist=nan(length(path3.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path3.dat)-1);path3.dist(i) = deg2km(distance(path3.dat(i,1),path3.dat(i,2),path3.dat(i+1,1),path3.dat(i+1,2)));end;
path3.speed=path3.dist./path3.dif;%figure,histogram(path1.speed,100)

outpath3=[repmat(3,length(pth3(:,1)),1),pth3(:,4:6),[0;path3.dist],[0;path3.speed]];

% PATH 4
path4.dat=sortrows(pth4(:,4:6),3,'descend');
path4.dif=abs(diff(path4.dat(:,3)));% time difference beteeen grid cell
path4.cum=sum(path4.dat(:,3));%total time along the trajectory
path4.dist=nan(length(path4.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path4.dat)-1);path4.dist(i) = deg2km(distance(path4.dat(i,1),path4.dat(i,2),path4.dat(i+1,1),path4.dat(i+1,2)));end;
path4.speed=path4.dist./path4.dif;%figure,histogram(path4.speed,100)

outpath4=[repmat(4,length(pth4(:,1)),1),pth4(:,4:6),[0;path4.dist],[0;path4.speed]];

% PATH 5
path5.dat=sortrows(pth5(:,4:6),3,'descend');
path5.dif=abs(diff(path5.dat(:,3)));% time difference beteeen grid cell
path5.cum=sum(path5.dat(:,3));%total time along the trajectory
path5.dist=nan(length(path5.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path5.dat)-1);path5.dist(i) = deg2km(distance(path5.dat(i,1),path5.dat(i,2),path5.dat(i+1,1),path5.dat(i+1,2)));end;
path5.speed=path5.dist./path5.dif;%figure,histogram(path5.speed,100)

outpath5=[repmat(5,length(pth5(:,1)),1),pth5(:,4:6),[0;path5.dist],[0;path5.speed]];

% PATH 6
path6.dat=sortrows(pth6(:,4:6),3,'descend');
path6.dif=abs(diff(path6.dat(:,3)));% time difference beteeen grid cell
path6.cum=sum(path6.dat(:,3));%total time along the trajectory
path6.dist=nan(length(path6.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path6.dat)-1);path6.dist(i) = deg2km(distance(path6.dat(i,1),path6.dat(i,2),path6.dat(i+1,1),path6.dat(i+1,2)));end;
path6.speed=path6.dist./path6.dif;%figure,histogram(path6.speed,100)

outpath6=[repmat(6,length(pth6(:,1)),1),pth6(:,4:6),[0;path6.dist],[0;path6.speed]];

out=[outpath1;outpath2;outpath3;outpath4;outpath5;outpath6];nout=length(out(:,1));
gen=load('Humancolonisation_FSTdata_v2.txt');out=[out,nan(nout,1)];%load genetic data
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg')
%coast = load('coastlines');seaclr1=hex2rgb('F0F8FF');%convert hexadecimal color into RGB triplet
%figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
%                    geoshow(coast.coastlat,coast.coastlon,'DisplayType','polygon','Facecolor','w'),
%                    geoshow(flipud(coast.coastlat),flipud(coast.coastlon),'DisplayType','polygon','FaceColor',seaclr1),
%                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
%                    geoshow(land, 'FaceColor', [0.8 0.8 0.8]);hold on,alpha(0.7);
%                    hold on,surfm(yg+0.5,xg2-0.5,gen);
for i=1:nout;
    vlat=abs(yg-out(i,2));%idlat=find(vlat==min(vlat));
    vlon=abs(xg2-out(i,3));%idlon=find(vlon==min(vlon));
    out(i,end)=gen(vlat==min(vlat),vlon==min(vlon));
    clear vlat vlon;% idlat idlon;
end;

% Now create the table
FigSD = table(out(:,1), out(:,2), out(:,3), out(:,4), out(:,5), out(:,6), out(:,7),'VariableNames', {'Trajectory', 'Latitude', 'Longitude', 'Date', 'Length', 'Speed', 'FST'});
% Save the table to a CSV file
writetable(FigSD, 'SourceData(AfricaPortugal).csv');

%saving the outcome
AfrPor.path1=path1;AfrPor.path2=path2;AfrPor.path3=path3;AfrPor.path4=path4;AfrPor.path5=path5;AfrPor.path6=path6;
save('AfricaPortugal.mat', '-struct', 'AfrPor');



%%=========================================================================
%% MEXICO-BRAZIL
%% the two best trajectoris are here #6 and #1
clear all, close all, clc,
pth1=load('Trajectory_Mexico-Brazil_(path#1).txt');pth2=load('Trajectory_Mexico-Brazil_(path#2).txt');pth3=load('Trajectory_Mexico-Brazil_(path#3).txt');pth4=load('Trajectory_Mexico-Brazil_(path#4).txt');
pth5=load('Trajectory_Mexico-Brazil_(path#5).txt');pth6=load('Trajectory_Mexico-Brazil_(path#6).txt');

% PATH 1
path1.dat=sortrows(pth1(:,4:6),3,'descend');
path1.dif=abs(diff(path1.dat(:,3)));% time difference beteeen grid cell
path1.cum=sum(path1.dat(:,3));%total time along the trajectory
path1.dist=nan(length(path1.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path1.dat)-1);path1.dist(i) = deg2km(distance(path1.dat(i,1),path1.dat(i,2),path1.dat(i+1,1),path1.dat(i+1,2)));end;
path1.speed=path1.dist./path1.dif;%figure,histogram(path1.speed,100)

outpath1=[ones(length(pth1(:,1)),1),pth1(:,4:6),[0;path1.dist],[0;path1.speed]];

% PATH 2
path2.dat=sortrows(pth2(:,4:6),3,'descend');
path2.dif=abs(diff(path2.dat(:,3)));% time difference beteeen grid cell
path2.cum=sum(path2.dat(:,3));%total time along the trajectory
path2.dist=nan(length(path2.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path2.dat)-1);path2.dist(i) = deg2km(distance(path2.dat(i,1),path2.dat(i,2),path2.dat(i+1,1),path2.dat(i+1,2)));end;
path2.speed=path2.dist./path2.dif;%figure,histogram(path1.speed,100)

outpath2=[repmat(2,length(pth2(:,1)),1),pth2(:,4:6),[0;path2.dist],[0;path2.speed]];

% PATH 3
path3.dat=sortrows(pth3(:,4:6),3,'descend');
path3.dif=abs(diff(path3.dat(:,3)));% time difference beteeen grid cell
path3.cum=sum(path3.dat(:,3));%total time along the trajectory
path3.dist=nan(length(path3.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path3.dat)-1);path3.dist(i) = deg2km(distance(path3.dat(i,1),path3.dat(i,2),path3.dat(i+1,1),path3.dat(i+1,2)));end;
path3.speed=path3.dist./path3.dif;%figure,histogram(path1.speed,100)

outpath3=[repmat(3,length(pth3(:,1)),1),pth3(:,4:6),[0;path3.dist],[0;path3.speed]];

% PATH 4
path4.dat=sortrows(pth4(:,4:6),3,'descend');
path4.dif=abs(diff(path4.dat(:,3)));% time difference beteeen grid cell
path4.cum=sum(path4.dat(:,3));%total time along the trajectory
path4.dist=nan(length(path4.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path4.dat)-1);path4.dist(i) = deg2km(distance(path4.dat(i,1),path4.dat(i,2),path4.dat(i+1,1),path4.dat(i+1,2)));end;
path4.speed=path4.dist./path4.dif;%figure,histogram(path4.speed,100)

outpath4=[repmat(4,length(pth4(:,1)),1),pth4(:,4:6),[0;path4.dist],[0;path4.speed]];

% PATH 5
path5.dat=sortrows(pth5(:,4:6),3,'descend');
path5.dif=abs(diff(path5.dat(:,3)));% time difference beteeen grid cell
path5.cum=sum(path5.dat(:,3));%total time along the trajectory
path5.dist=nan(length(path5.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path5.dat)-1);path5.dist(i) = deg2km(distance(path5.dat(i,1),path5.dat(i,2),path5.dat(i+1,1),path5.dat(i+1,2)));end;
path5.speed=path5.dist./path5.dif;%figure,histogram(path5.speed,100)

outpath5=[repmat(5,length(pth5(:,1)),1),pth5(:,4:6),[0;path5.dist],[0;path5.speed]];

% PATH 6
path6.dat=sortrows(pth6(:,4:6),3,'descend');
path6.dif=abs(diff(path6.dat(:,3)));% time difference beteeen grid cell
path6.cum=sum(path6.dat(:,3));%total time along the trajectory
path6.dist=nan(length(path6.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path6.dat)-1);path6.dist(i) = deg2km(distance(path6.dat(i,1),path6.dat(i,2),path6.dat(i+1,1),path6.dat(i+1,2)));end;
path6.speed=path6.dist./path6.dif;%figure,histogram(path6.speed,100)

outpath6=[repmat(6,length(pth6(:,1)),1),pth6(:,4:6),[0;path6.dist],[0;path6.speed]];

out=[outpath1;outpath2;outpath3;outpath4;outpath5;outpath6];nout=length(out(:,1));
gen=load('Humancolonisation_FSTdata_v2.txt');out=[out,nan(nout,1)];%load genetic data
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg')
%coast = load('coastlines');seaclr1=hex2rgb('F0F8FF');%convert hexadecimal color into RGB triplet
%figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
%                    geoshow(coast.coastlat,coast.coastlon,'DisplayType','polygon','Facecolor','w'),
%                    geoshow(flipud(coast.coastlat),flipud(coast.coastlon),'DisplayType','polygon','FaceColor',seaclr1),
%                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
%                    geoshow(land, 'FaceColor', [0.8 0.8 0.8]);hold on,alpha(0.7);
%                    hold on,surfm(yg+0.5,xg2-0.5,gen);
for i=1:nout;
    vlat=abs(yg-out(i,2));%idlat=find(vlat==min(vlat));
    vlon=abs(xg2-out(i,3));%idlon=find(vlon==min(vlon));
    out(i,end)=gen(vlat==min(vlat),vlon==min(vlon));
    clear vlat vlon;% idlat idlon;
end;

% Now create the table
FigSD = table(out(:,1), out(:,2), out(:,3), out(:,4), out(:,5), out(:,6), out(:,7),'VariableNames', {'Trajectory', 'Latitude', 'Longitude', 'Date', 'Length', 'Speed', 'FST'});
% Save the table to a CSV file
writetable(FigSD, 'SourceData(MexicoBrazil).csv');

%saving the outcome
MexBra.path1=path1;MexBra.path2=path2;MexBra.path3=path3;MexBra.path4=path4;MexBra.path5=path5;MexBra.path6=path6;
save('MexicoBrazil.mat', '-struct', 'MexBra');


%%=========================================================================
%% MEXICO-BERINGIA
%% the two best trajectoris are here #6 and #1
clear all, close all, clc,   
pth1=load('Trajectory_Bering-Mexico_(path#1).txt');pth2=load('Trajectory_Bering-Mexico_(path#2).txt');pth3=load('Trajectory_Bering-Mexico_(path#3).txt');
pth4=load('Trajectory_Bering-Mexico_(path#4).txt');pth5=load('Trajectory_Bering-Mexico_(path#5).txt');pth6=load('Trajectory_Bering-Mexico_(path#6).txt');
pth7=load('Trajectory_Bering-Mexico_(path#7).txt');pth8=load('Trajectory_Bering-Mexico_(path#8).txt');

% PATH 1
path1.dat=sortrows(pth1(:,4:6),3,'descend');
path1.dif=abs(diff(path1.dat(:,3)));% time difference beteeen grid cell
path1.cum=sum(path1.dat(:,3));%total time along the trajectory
path1.dist=nan(length(path1.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path1.dat)-1);path1.dist(i) = deg2km(distance(path1.dat(i,1),path1.dat(i,2),path1.dat(i+1,1),path1.dat(i+1,2)));end;
path1.speed=path1.dist./path1.dif;%figure,histogram(path1.speed,100)

outpath1=[ones(length(pth1(:,1)),1),pth1(:,4:6),[0;path1.dist],[0;path1.speed]];

% PATH 2
path2.dat=sortrows(pth2(:,4:6),3,'descend');
path2.dif=abs(diff(path2.dat(:,3)));% time difference beteeen grid cell
path2.cum=sum(path2.dat(:,3));%total time along the trajectory
path2.dist=nan(length(path2.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path2.dat)-1);path2.dist(i) = deg2km(distance(path2.dat(i,1),path2.dat(i,2),path2.dat(i+1,1),path2.dat(i+1,2)));end;
path2.speed=path2.dist./path2.dif;%figure,histogram(path1.speed,100)

outpath2=[repmat(2,length(pth2(:,1)),1),pth2(:,4:6),[0;path2.dist],[0;path2.speed]];

% PATH 3
path3.dat=sortrows(pth3(:,4:6),3,'descend');
path3.dif=abs(diff(path3.dat(:,3)));% time difference beteeen grid cell
path3.cum=sum(path3.dat(:,3));%total time along the trajectory
path3.dist=nan(length(path3.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path3.dat)-1);path3.dist(i) = deg2km(distance(path3.dat(i,1),path3.dat(i,2),path3.dat(i+1,1),path3.dat(i+1,2)));end;
path3.speed=path3.dist./path3.dif;%figure,histogram(path1.speed,100)

outpath3=[repmat(3,length(pth3(:,1)),1),pth3(:,4:6),[0;path3.dist],[0;path3.speed]];

% PATH 4
path4.dat=sortrows(pth4(:,4:6),3,'descend');
path4.dif=abs(diff(path4.dat(:,3)));% time difference beteeen grid cell
path4.cum=sum(path4.dat(:,3));%total time along the trajectory
path4.dist=nan(length(path4.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path4.dat)-1);path4.dist(i) = deg2km(distance(path4.dat(i,1),path4.dat(i,2),path4.dat(i+1,1),path4.dat(i+1,2)));end;
path4.speed=path4.dist./path4.dif;%figure,histogram(path4.speed,100)

outpath4=[repmat(4,length(pth4(:,1)),1),pth4(:,4:6),[0;path4.dist],[0;path4.speed]];

% PATH 5
path5.dat=sortrows(pth5(:,4:6),3,'descend');
path5.dif=abs(diff(path5.dat(:,3)));% time difference beteeen grid cell
path5.cum=sum(path5.dat(:,3));%total time along the trajectory
path5.dist=nan(length(path5.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path5.dat)-1);path5.dist(i) = deg2km(distance(path5.dat(i,1),path5.dat(i,2),path5.dat(i+1,1),path5.dat(i+1,2)));end;
path5.speed=path5.dist./path5.dif;%figure,histogram(path5.speed,100)

outpath5=[repmat(5,length(pth5(:,1)),1),pth5(:,4:6),[0;path5.dist],[0;path5.speed]];

% PATH 6
path6.dat=sortrows(pth6(:,4:6),3,'descend');
path6.dif=abs(diff(path6.dat(:,3)));% time difference beteeen grid cell
path6.cum=sum(path6.dat(:,3));%total time along the trajectory
path6.dist=nan(length(path6.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path6.dat)-1);path6.dist(i) = deg2km(distance(path6.dat(i,1),path6.dat(i,2),path6.dat(i+1,1),path6.dat(i+1,2)));end;
path6.speed=path6.dist./path6.dif;%figure,histogram(path6.speed,100)

outpath6=[repmat(6,length(pth6(:,1)),1),pth6(:,4:6),[0;path6.dist],[0;path6.speed]];

% PATH 7
path7.dat=sortrows(pth7(:,4:6),3,'descend');
path7.dif=abs(diff(path7.dat(:,3)));% time difference beteeen grid cell
path7.cum=sum(path7.dat(:,3));%total time along the trajectory
path7.dist=nan(length(path7.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path7.dat)-1);path7.dist(i) = deg2km(distance(path7.dat(i,1),path7.dat(i,2),path7.dat(i+1,1),path7.dat(i+1,2)));end;
path7.speed=path7.dist./path7.dif;%figure,histogram(path7.speed,100)

outpath7=[repmat(7,length(pth7(:,1)),1),pth7(:,4:6),[0;path7.dist],[0;path7.speed]];

% PATH 8
path8.dat=sortrows(pth8(:,4:6),3,'descend');
path8.dif=abs(diff(path8.dat(:,3)));% time difference beteeen grid cell
path8.cum=sum(path8.dat(:,3));%total time along the trajectory
path8.dist=nan(length(path8.dat)-1,1);%calculation distance between gridcell along the trajectories
for i=1:(length(path8.dat)-1);path8.dist(i) = deg2km(distance(path8.dat(i,1),path8.dat(i,2),path8.dat(i+1,1),path8.dat(i+1,2)));end;
path8.speed=path8.dist./path8.dif;%figure,histogram(path8.speed,100)

outpath8=[repmat(8,length(pth8(:,1)),1),pth8(:,4:6),[0;path8.dist],[0;path8.speed]];

out=[outpath1;outpath2;outpath3;outpath4;outpath5;outpath6;outpath7;outpath8];nout=length(out(:,1));
gen=load('Humancolonisation_FSTdata_v2.txt');out=[out,nan(nout,1)];%load genetic data
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg')
%coast = load('coastlines');seaclr1=hex2rgb('F0F8FF');%convert hexadecimal color into RGB triplet
%figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
%                    geoshow(coast.coastlat,coast.coastlon,'DisplayType','polygon','Facecolor','w'),
%                    geoshow(flipud(coast.coastlat),flipud(coast.coastlon),'DisplayType','polygon','FaceColor',seaclr1),
%                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
%                    geoshow(land, 'FaceColor', [0.8 0.8 0.8]);hold on,alpha(0.7);
%                    hold on,surfm(yg+0.5,xg2-0.5,gen);
for i=1:nout;
    vlat=abs(yg-out(i,2));%idlat=find(vlat==min(vlat));
    vlon=abs(xg2-out(i,3));%idlon=find(vlon==min(vlon));
    out(i,end)=gen(vlat==min(vlat),vlon==min(vlon));
    clear vlat vlon;% idlat idlon;
end;

% Now create the table
FigSD = table(out(:,1), out(:,2), out(:,3), out(:,4), out(:,5), out(:,6), out(:,7),'VariableNames', {'Trajectory', 'Latitude', 'Longitude', 'Date', 'Length', 'Speed', 'FST'});
% Save the table to a CSV file
writetable(FigSD, 'SourceData(BeringiaMexico).csv');


%saving the outcome
BerMex.path1=path1;BerMex.path2=path2;BerMex.path3=path3;BerMex.path4=path4;BerMex.path5=path5;BerMex.path6=path6;BerMex.path7=path7;BerMex.path8=path8;
save('BeringiaMexico.mat', '-struct', 'BerMex');


                    