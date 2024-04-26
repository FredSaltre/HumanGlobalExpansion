program Combined_maps


%%=========================================================================
%% STEP #0: PLOT OF THE HUMAN ARRIVAL ONLY
%%=========================================================================
clear all, close all, clc,
mat=load('maskedHumantimingraster_v2.txt');
matsd=load('maskedHumantimingraster(StandardDeviation)_v2.txt');

%%=========================================================================
%% STEP #1: REPLOT OF THE ANOMALIE OF TEMPERATURE : exit - arrival
%%=========================================================================
% temperature
clear all, close all, clc,
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg')
base=load('GlobalClimateAnomalie(pre80ka)_temperature.txt');base2=[base(:,tr:end),base(:,1:tr-1)];
pred=load('NewTempAtArrivalTiming.txt');outT=pred-base2;outT(1:48,315:end)=NaN;
icesheet1=load('icesheet16ka.txt');icesheet2=load('IceSheetAntarticaGreenlandEurope.txt');
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);Cice=brewermap(13,'Greys');
clr2=brewermap(13,'Reds');clr3=brewermap(13,'Blues');clr4=[flipud(clr2);clr3];clr4(1:4,:)=[];clr4(end-4:end,:)=[];
rsys=load('worldmainrivers.txt');vrsys=unique(rsys(:,1));nr=length(vrsys);

figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
                    geoshow(land, 'FaceColor', [0.8 0.8 0.8])
                    hold on,surfm(yg+0.5,xg2-0.5,outT);colormap(flipud(clr4));caxis([-3,3]);%hold on,alpha(0.6)
                    hold off,geoshow(icesheet1(:,2),icesheet1(:,1),'DisplayType','polygon','Facecolor','w','Edgecolor','w'); 
                    hold on,plotm(icesheet2(:,2),icesheet2(:,1),'ow','Marker','square','Markersize',2.2,'MarkerFacecolor','w','MarkerEdgecolor','w')
                    for i=1:nr;id=find(rsys(:,1)==vrsys(i));hold on,plotm(rsys(id,3),rsys(id,2),'-b','color',seaclr(15,:),'Linewidth',0.5);clear id;end;
                    hold on,plotm(flipud(coast.lat),flipud(coast.long),'-k');
                    tr=title('Temperature changes at the timing of human arrival (C)');set(tr,'Fontsize',13);
save -ascii Final_AnomalieTemperature(exit-Arrival).txt outT;

% precipitation
clear all, close all, clc,
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg')
base=load('GlobalClimateAnomalie(pre80ka)_precipitation.txt');base2=[base(:,tr:end),base(:,1:tr-1)];
pred=load('NewPrecipAtArrivalTiming.txt');icesheet1=load('icesheet16ka.txt');icesheet2=load('IceSheetAntarticaGreenlandEurope.txt');
outP=pred-base2;outP(1:48,315:end)=NaN;

coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);Cice=brewermap(13,'Greys');
clr2=brewermap(13,'YlGn');%clr3=brewermap(13,'Blues');clr4=[flipud(clr2);clr3];clr4(1:4,:)=[];clr4(end-4:end,:)=[];
rsys=load('worldmainrivers.txt');vrsys=unique(rsys(:,1));nr=length(vrsys);


figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
                    geoshow(land, 'FaceColor', [0.8 0.8 0.8])
                    hold on,surfm(yg+0.5,xg2-0.5,outP);colormap(clr2);caxis([-0.5,0.5]);%hold on,alpha(0.6)
                    hold off,geoshow(icesheet1(:,2),icesheet1(:,1),'DisplayType','polygon','Facecolor','w','Edgecolor','w'); 
                    hold on,plotm(icesheet2(:,2),icesheet2(:,1),'ow','Marker','square','Markersize',2.2,'MarkerFacecolor','w','MarkerEdgecolor','w')
                    for i=1:nr;id=find(rsys(:,1)==vrsys(i));hold on,plotm(rsys(id,3),rsys(id,2),'-b','color',seaclr(15,:),'Linewidth',0.5);clear id;end;
                    hold on,plotm(flipud(coast.lat),flipud(coast.long),'-k');
                    tr=title('Precipitation changes at the timing of human arrival (C)');set(tr,'Fontsize',13);
save -ascii Final_AnomaliePrecipitation(exit-Arrival).txt outP;
                    
% vegetation
% if from forest (1) to non forest (2) => 20 => light orange
%         non forest (2) to forest (1) => 5 => light green
%         forest (1) to forest (1) => 15 => dark green
%         non forest (2) to non forest (2) = 10 => dark orange
clear all, close all, clc,
base=load('GlobalVegetation(pre80ka).txt');pred=load('NewForestNonForestArrivalTiming.txt');[r,l]=size(base);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg')
base2=[base(:,tr:end),base(:,1:tr-1)];TF=isnan(pred);[row,col]=find(TF==0);out=nan(r,l);nr=length(row);h = waitbar(0,'Please wait...');
for i=1:nr;waitbar(i/nr)
    if (base2(row(i),col(i))==1 & pred(row(i),col(i))==2);out(row(i),col(i))=20;
    elseif (base2(row(i),col(i))==2 & pred(row(i),col(i))==1);out(row(i),col(i))=5;
    elseif (base2(row(i),col(i))==1 & pred(row(i),col(i))==1);out(row(i),col(i))=15;
    else (base2(row(i),col(i))==2 & pred(row(i),col(i))==2);out(row(i),col(i))=10; 
    end;
end;close(h);out(1:48,315:end)=NaN;
icesheet1=load('icesheet16ka.txt');icesheet2=load('IceSheetAntarticaGreenlandEurope.txt');
rsys=load('worldmainrivers.txt');vrsys=unique(rsys(:,1));nr=length(vrsys);

coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);Cice=brewermap(13,'Greys');
clr2=brewermap(12,'Paired');clr3=[clr2(3,:);clr2(8,:);clr2(4,:);clr2(7,:)];


figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
                    geoshow(land, 'FaceColor', [0.8 0.8 0.8])
                    hold on,surfm(yg+0.5,xg2-0.5,out);colormap(clr3);hold on,alpha(0.5)
                    hold off,geoshow(icesheet1(:,2),icesheet1(:,1),'DisplayType','polygon','Facecolor','w','Edgecolor','w'); 
                    hold on,plotm(icesheet2(:,2),icesheet2(:,1),'ow','Marker','square','Markersize',2.2,'MarkerFacecolor','w','MarkerEdgecolor','w')
                    for i=1:nr;id=find(rsys(:,1)==vrsys(i));hold on,plotm(rsys(id,3),rsys(id,2),'-b','color',seaclr(15,:),'Linewidth',0.5);clear id;end;
                    hold on,plotm(flipud(coast.lat),flipud(coast.long),'-k');
                    tr=title('Vegetation changes at the timing of human arrival (C)');set(tr,'Fontsize',13);                   
save -ascii Final_Forestchanges(exit-Arrival).txt out;

%%=========================================================================
%% STEP #2: COMBINING THE DIFFERENT POSSIBILITY OF CHANGES 
%%=========================================================================
clear all, close all, clc,
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');
T=load('Final_AnomalieTemperature(exit-Arrival).txt');P=load('Final_AnomaliePrecipitation(exit-Arrival).txt');
V=load('Final_Forestchanges(exit-Arrival).txt');DR=load('Distance2nearestRiver(raster).txt');[l,c]=size(T);
TF=isnan(T);[row,col]=find(TF==0);out=nan(l,c);nr=length(row);h = waitbar(0,'Please wait...');
for i=1:nr;waitbar(i/nr)
    if (T(row(i),col(i)) >= 0); %temperature increased
        if (P(row(i),col(i)) >= 0)%precipitation increased
            switch V(row(i),col(i))
                case 5%vegetation = forest
                    out(row(i),col(i))=1;
                case 15
                    out(row(i),col(i))=1;
                otherwise %vegetation = grassland
                    out(row(i),col(i))=3;
            end;
        else%precipitation decreased
           switch V(row(i),col(i))
                case 5%vegetation = forest
                    out(row(i),col(i))=2;
                case 15%vegetation = forest
                    out(row(i),col(i))=2;
                otherwise%vegetation = grassland
                    out(row(i),col(i))=4;
            end; 
        end
    else %temperature decreased
        if (P(row(i),col(i)) >= 0)%precipitation increased
            switch V(row(i),col(i))
                case 5%vegetation = forest
                    out(row(i),col(i))=5;
                case 15
                    out(row(i),col(i))=5;
                otherwise %vegetation = grassland
                    out(row(i),col(i))=7;
            end;
        else%precipitation decreased
           switch V(row(i),col(i))
                case 5%vegetation = forest
                    out(row(i),col(i))=6;
                case 15%vegetation = forest
                    out(row(i),col(i))=6;
                otherwise%vegetation = grassland
                    out(row(i),col(i))=8;
            end; 
        end                
    end;
end;close(h);
res=[];for i=1:c;res=[res;out(:,i)];end;TF=isnan(res);id=find(TF==1);res(id)=[];
DR2=(DR-min(min(DR)))/(max(max(DR))-min(min(DR)));
out2=nan(l,c);h = waitbar(0,'Please wait...');
for i=1:nr;waitbar(i/nr);out2(row(i),col(i))=(out(row(i),col(i))-1)+DR2(row(i),col(i));end;close(h);

%%% NEED HERE TO CONVERT THE RASTER VALUE IN LON LAT + T, P, VEG, DIST RIV, FINAL CODE 
SDat=nan(nr,8);%prepare the Source Data file output [,1] = lat, [,2] = lon, [,3] = T, [,4] = P, [,5] = V, [,6]= Dist river, [,7]= rescale Dist river, [,8] = colour code for figure 1
% vegetation colour code
% if from forest (1) to non forest (2) => 20 => light orange
%         non forest (2) to forest (1) => 5 => light green
%         forest (1) to forest (1) => 15 => dark green
%         non forest (2) to non forest (2) = 10 => dark orange
h = waitbar(0,'Please wait...');
for i=1:nr;waitbar(i/nr);
    SDat(i,:)=[yg(row(i)),xg2(col(i)),T(row(i),col(i)),P(row(i),col(i)),V(row(i),col(i)),DR(row(i),col(i)),DR2(row(i),col(i)),out2(row(i),col(i))];
end;close(h);


coast = load('coastlines');seaclr1=hex2rgb('F0F8FF');%convert hexadecimal color into RGB triplet
clr1=brewermap(13,'Blues');clr2=brewermap(13,'Reds');clr3=brewermap(13,'Greens');
t1=brewermap(13,'Oranges');t2=brewermap(13,'YlOrBr');
clr6=[t2(1:3,:);t1(1:8,:)];
clr=[clr1(3:13,:);clr6;clr2(3:13,:);clr3(3:13,:)];

figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.coastlat,coast.coastlon,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.coastlat),flipud(coast.coastlon),'DisplayType','polygon','FaceColor',seaclr1),
                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
                    geoshow(land, 'FaceColor', [0.8 0.8 0.8]);hold on,alpha(0.7);
                    hold on,scatterm(SDat(:,1)+0.5,SDat(:,2)-0.5,5,SDat(:,end));colormap(clr);clim([4 8]);%

Latitude = SDat(:,3);Longitude = SDat(:,2); T_anomaly = SDat(:,3); P_anomaly = SDat(:,4); Veg_ID = SDat(:,5); Dist2Riv = SDat(:,6); Scale_Dist2Riv = SDat(:,7); Fig_code = SDat(:,8);
% Now create the table
Fig1_SD = table(Latitude,Longitude, T_anomaly, P_anomaly, Veg_ID, Dist2Riv, Scale_Dist2Riv, Fig_code,'VariableNames', {'Latitude', 'Longitude', 'T_anomaly', 'P_anomaly', 'Veg_ID', 'Dist2Riv', 'Scale_Dist2Riv', 'Fig_code'});

% Save the table to a CSV file
writetable(Fig1_SD, 'SourceData_Fig1.csv');




icesheet1=load('icesheet16ka.txt');icesheet2=load('IceSheetAntarticaGreenlandEurope.txt');

coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);Cice=brewermap(13,'Greys');
rsys=load('worldmainrivers.txt');vrsys=unique(rsys(:,1));nr2=length(vrsys);


figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
                    geoshow(land, 'FaceColor', [0.8 0.8 0.8]);hold on,alpha(0.7);
                    hold on,surfm(yg+0.5,xg2-0.5,out2);colormap(clr);caxis([4 8]);%
                    hold off,geoshow(icesheet1(:,2),icesheet1(:,1),'DisplayType','polygon','Facecolor','w','Edgecolor','w'); 
                    hold on,plotm(icesheet2(:,2),icesheet2(:,1),'ow','Marker','square','Markersize',2.2,'MarkerFacecolor','w','MarkerEdgecolor','w')
                    for i=1:nr2;id=find(rsys(:,1)==vrsys(i));hold on,plotm(rsys(id,3),rsys(id,2),'-b','color',seaclr(15,:),'Linewidth',0.5);clear id;end;
                    hold on,plotm(flipud(coast.lat),flipud(coast.long),'-k');
                    tr=title('Vegetation changes at the timing of human arrival (C)');set(tr,'Fontsize',13);       



