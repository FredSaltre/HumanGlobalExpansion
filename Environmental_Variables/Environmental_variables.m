program Environmental_variables

%%=========================================================================
%% RETURN A RASTER WITH THE DISTANCE OF EACH GRIDCELL TO THE NEAREST COAST
%%=========================================================================
clear all, close all, clc,
mat=load('maskedHumantimingraster_v2.txt');[r,c]=size(mat);tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');
TF=isnan(mat);[row,col]=find(TF==0);nr=length(row);out=[row,col,nan(nr,3)];for i=1:nr;out(i,3:4)=[yg(row(i)),xg2(col(i))];end;
[row2,col2]=find(TF==1);nr2=length(row2);coastal=[row2,col2,nan(nr2,2)];for i=1:nr2;coastal(i,3:4)=[yg(row2(i)),xg2(col2(i))];end;
h = waitbar(0,'Please wait...');
for i=1:nr;waitbar(i/nr)
    [vc,~]=distance(out(i,3),out(i,4),coastal(:,3),coastal(:,4));
    id=find(vc==min(vc));out(i,5)=deg2km(vc(id(1)));
    clear vc id;
end;close(h);
res=nan(r,c);for i=1:nr;res(out(i,1),out(i,2))=out(i,5);end;
save -ascii NewDistance2coast.txt res;

icesheet1=load('icesheet16ka.txt');icesheet2=load('IceSheetAntarticaGreenlandEurope.txt');
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);clr3=brewermap(13,'YlOrBr');%Cice=brewermap(13,'Greys');
figure('Color','w'),worldmap({'world'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
                    geoshow(land, 'FaceColor', [0.8 0.8 0.8]);
                    %hold on;plotm(riv(:,3),riv(:,2),'.b'); %hold on;plotm(chum(:,3),chum(:,4),'.r');
                    
                    hold on,surfm(yg+0.5,xg2-0.5,res);colormap(clr3);
                    %hold off,geoshow(icesheet1(:,2),icesheet1(:,1),'DisplayType','polygon','Facecolor',Cice(7,:),'Edgecolor',Cice(7,:)); 
                    hold off,geoshow(icesheet1(:,2),icesheet1(:,1),'DisplayType','polygon','Facecolor','w','Edgecolor','w'); 
                    %hold on,geoshow(icesheet2(:,2),icesheet2(:,1),'DisplayType','point','Marker','square','Markersize',2,'MarkerFacecolor',Cice(7,:),'MarkerEdgecolor',Cice(7,:)),
                    hold on,geoshow(icesheet2(:,2),icesheet2(:,1),'DisplayType','point','Marker','square','Markersize',2.2,'MarkerFacecolor','w','MarkerEdgecolor','w'),
                    hold on,plotm(flipud(coast.lat),flipud(coast.long),'-k');tr=title('Distance to nearest coast');set(tr,'Fontsize',14);

%%=========================================================================
%% RETURN A RASTER WITH THE DISTANCE OF EACH GRIDCELL TO THE NEAREST RIVER 
%%=========================================================================
clear all, close all, clc,
rsys=load('worldmainrivers.txt');vrsys=unique(rsys(:,1));
hum=load('maskedHumantimingraster_v2.txt');TF=isnan(hum);[row,col]=find(TF==0);chum=nan(length(row),4);
for i=1:length(row);chum(i,:)=[row(i),col(i),yg(row(i)),xg2(col(i))];end;nc=length(chum(:,1));chum=[chum,nan(nc,1)];
for i=1:nc;
    [vc,~]=distance(chum(i,3),chum(i,4),riv(:,3),riv(:,2));
    id=find(vc==min(vc));chum(i,5)=deg2km(vc(id(1)));
    clear vc id;
end;
out=nan(180,360);for i=1:nc;out(chum(i,1),chum(i,2))=chum(i,5);end;
save -ascii Distance2nearestRiver(raster).txt out;

out=load('Distance2nearestRiver(raster).txt');
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');
icesheet1=load('icesheet16ka.txt');icesheet2=load('IceSheetAntarticaGreenlandEurope.txt');
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);clr3=flipud(brewermap(13,'Blues'));clr2=clr3(1,:);clr3(1:6,:)=[];%Cice=brewermap(13,'Greys');


figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
                    geoshow(land, 'FaceColor', [0.8 0.8 0.8])
                    hold on,surfm(yg+0.5,xg2-0.5,out);colormap(clr3);alpha 0.8;%caxis([-0 2])
                    hold off,geoshow(icesheet1(:,2),icesheet1(:,1),'DisplayType','polygon','Facecolor','w','Edgecolor','w'); 
                    hold on,plotm(icesheet2(:,2),icesheet2(:,1),'ow','Marker','square','Markersize',2.2,'MarkerFacecolor','w','MarkerEdgecolor','w')
                    for i=1:nr;id=find(rsys(:,1)==vrsys(i));hold on,plotm(rsys(id,3),rsys(id,2),'-b','color',clr2,'Linewidth',0.5);clear id;end;
                    hold on,plotm(flipud(coast.lat),flipud(coast.long),'-k');
                    tr=title('Distance to the nearest river');set(tr,'Fontsize',14);


figure('Color','w'),worldmap({'world'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
                    geoshow(land, 'FaceColor', [0.8 0.8 0.8]);
                    %hold on;plotm(riv(:,3),riv(:,2),'.b'); %hold on;plotm(chum(:,3),chum(:,4),'.r');
                    
                    hold on,surfm(yg+0.5,xg2-0.5,out);colormap(clr3);alpha 0.8;%caxis([-0 2])
                    %hold off,geoshow(icesheet1(:,2),icesheet1(:,1),'DisplayType','polygon','Facecolor',Cice(7,:),'Edgecolor',Cice(7,:)); 
                    hold off,geoshow(icesheet1(:,2),icesheet1(:,1),'DisplayType','polygon','Facecolor','w','Edgecolor','w'); 
                    %hold on,geoshow(icesheet2(:,2),icesheet2(:,1),'DisplayType','point','Marker','square','Markersize',2,'MarkerFacecolor',Cice(7,:),'MarkerEdgecolor',Cice(7,:)),
                    hold on,geoshow(icesheet2(:,2),icesheet2(:,1),'DisplayType','point','Marker','square','Markersize',2.2,'MarkerFacecolor','w','MarkerEdgecolor','w'),
                    for i=1:nr;id=find(rsys(:,1)==vrsys(i));hold on,plotm(rsys(id,3),rsys(id,2),'-b','color',clr2,'Linewidth',0.5);clear id;end;
                    hold on,plotm(flipud(coast.lat),flipud(coast.long),'-k');tr=title('Distance to the nearest river');set(tr,'Fontsize',14);

%%=========================================================================
%% RETURN A RASTER WITH AN INDEX OF RUGGEDNESS FOR EACH GRID CELL WORLDWIDE 
%%=========================================================================
clear all, close all, clc,
mat.alt=ncread('elev.1-deg.nc','data');mat.lat=ncread('elev.1-deg.nc','lat');mat.lon=ncread('elev.1-deg.nc','lon');
mat.time=load('maskedHumantimingraster_v2.txt');tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');%fichier humain timing
mat.alt=mat.alt';mat.lon=[0.5:179.5,-179.5:1:-0.5];%fichier altitude -> a mettre sur la meme grille
id=find(mat.lon==xg2(1));mat.alt2=[mat.alt(:,id:end),mat.alt(:,1:id-1)];%on matche les deux grgille aux meme longitudes
nr=180;nc=360;cpt=0;mat.rug=nan(nr,nc);h = waitbar(0,'Please wait...');
for i=2:nr-1;
    for j=2:nc-1;ngh=nan(9,1);cpt=cpt+1;waitbar(cpt/(nr*nc))
        ngh(1,1)=(mat.alt2(i-1,j-1)-mat.alt2(i,j)).^2;
        ngh(2,1)=(mat.alt2(i-1,j)-mat.alt2(i,j)).^2;
        ngh(3,1)=(mat.alt2(i-1,j+1)-mat.alt2(i,j)).^2;
        ngh(4,1)=(mat.alt2(i,j-1)-mat.alt2(i,j)).^2;
        ngh(5,1)=(mat.alt2(i,j)-mat.alt2(i,j)).^2;
        ngh(6,1)=(mat.alt2(i,j+1)-mat.alt2(i,j)).^2;
        ngh(7,1)=(mat.alt2(i+1,j-1)-mat.alt2(i,j)).^2;
        ngh(8,1)=(mat.alt2(i+1,j)-mat.alt2(i,j)).^2;
        ngh(9,1)=(mat.alt2(i+1,j+1)-mat.alt2(i,j)).^2;
        mat.rug(i,j)=sqrt(sum(ngh));clear ngh;        
    end;
end;close(h);
TF=isnan(mat.time);[row,col]=find(TF==1);for i=1:length(row);mat.rug(row(i),col(i))=NaN;end;


out=mat.rug;save -ascii NewRuggedness.txt out;


