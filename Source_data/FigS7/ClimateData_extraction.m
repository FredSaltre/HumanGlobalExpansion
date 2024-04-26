program ClimateData_extraction

%%=========================================================================
%% STEP#1: DESIGN OF MASKS
%%=========================================================================
clear all, close all, clc,
mat1=load('Europe+Asie_raster.txt');Ctr=load('countrylim2.txt');
clon=-179.5:1:179.5;clat=-89.5:89.5;clat=flipud(clat');
%mask for Eastern Europe
idEE=[23;196;248;190;64;215;110;37;195;157];nEE=length(idEE);cIdEE=[];for i=1:nEE;[row,col]=find(Ctr==idEE(i));cIdEE=[cIdEE;[row,col]];clear row col;end;%cIdEE=[cIdEE(:,1)+6,cIdEE(:,2)-1];cIdEE(cIdEE(:,2)==0,2)=360;
nEE2=length(cIdEE(:,1));for i=1:nEE2;mat1(cIdEE(i,1),cIdEE(i,2))=6;end;
%mask for North America
idNA=[251,43];nNA=length(idNA);cIdNA=[];for i=1:nNA;[row,col]=find(Ctr==idNA(i));cIdNA=[cIdNA;[row,col]];clear row col;end;
nNA2=length(cIdNA(:,1));for i=1:nNA2;mat1(cIdNA(i,1),cIdNA(i,2))=8;end;
%mask for Central America
idCA=[155,61,69,105,100,171,72,59,184,245];nCA=length(idCA);cIdCA=[];for i=1:nCA;[row,col]=find(Ctr==idCA(i));cIdCA=[cIdCA;[row,col]];clear row col;end;
nCA2=length(cIdCA(:,1));for i=1:nCA2;mat1(cIdCA(i,1),cIdCA(i,2))=12;end;
%mask for South America
idSA=[53,255,104,226,84,34,70,187,29,48,10,186,252];nSA=length(idSA);cIdSA=[];for i=1:nSA;[row,col]=find(Ctr==idSA(i));cIdSA=[cIdSA;[row,col]];clear row col;end;
nSA2=length(cIdSA(:,1));for i=1:nSA2;mat1(cIdSA(i,1),cIdSA(i,2))=15;end;
%mask for Australia
idA=14;nA=length(idA);cIdA=[];for i=1:nA;[row,col]=find(Ctr==idA(i));cIdA=[cIdA;[row,col]];clear row col;end;
nA2=length(cIdA(:,1));for i=1:nA2;mat1(cIdA(i,1),cIdA(i,2))=18;end;

TF=isnan(mat1);[row,col]=find(TF==0);nr=length(row);cres=nan(nr,3);for i=1:nr;cres(i,:)=[row(i),col(i),mat1(row(i),col(i))];end;
cres(:,1:2)=[cres(:,1)+6,cres(:,2)-1];cres(cres(:,2)==0,2)=360;
res=nan(180,360);for i=1:nr;res(cres(i,1),cres(i,2))=cres(i,3);end;
save -ascii AllRegionRaster.txt res;

coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);clr=brewermap(12,'Set3');
figure('Color','w'),worldmap({'world'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(clat+0.5,clon+0.5,res);colormap(clr);caxis([5000,55000]);hold on,
                    hold on,plotm(flipud(coast.lat),flipud(coast.long),'-k')
                    
res2=load('AllRegionRaster_v2.txt');clr2=brewermap(8,'Set2');
figure('Color','w'),worldmap({'world'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(clat+0.5,clon+0.5,res2);colormap(clr2);hold on,
                    hold on,plotm(flipud(coast.lat),flipud(coast.long),'-k')                   
                    caxis([1,8]);
                    
%%=========================================================================
%% STEP#2: TEMPERATURES EXTRACTION BY REGIONS
%%=========================================================================
clear all, close all, clc,   
mat=load('AllRegionRaster_v2.txt');tps=[0,1000:1000:22000,24000:2000:80000,84000:4000:120000];nt=length(tps);
clon=-179.5:1:179.5;clat=-89.5:89.5;clat=flipud(clat');
TF=isnan(mat);[row,col]=find(TF==0);nr=length(row);cmat=nan(nr,nt+3);coord=nan(nr,2);
for i=1:nr;
    cmat(i,1:3)=[row(i),col(i),mat(row(i),col(i))];
    coord(i,:)=[clat(row(i)),clon(col(i))];
end;
h = waitbar(0,'Please wait...');
for t=1:nt;waitbar(t/nt)
    eval(['clim=load(''' num2str(tps(t)) '_amean_temp.txt'')']);
    for i=1:nr;cmat(i,t+3)=clim(cmat(i,1),cmat(i,2));end;
    clear clim;
end;close(h);save -ascii TempAnom_Allregion.txt cmat;


%%formatting for SOurce Data files
SDtp=load('TempAnom_Allregion.txt');
SDtpFigS7=[coord,SDtp(:,4:end)];

coast = load('coastlines');seaclr1=hex2rgb('F0F8FF');%convert hexadecimal color into RGB triplet

figure('Color','w'),ax=worldmap({'world'});setm(ax, 'Origin', [0 150 0]),gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.coastlat,coast.coastlon,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.coastlat),flipud(coast.coastlon),'DisplayType','polygon','FaceColor',seaclr1),
                    land = shaperead('landareas.shp', 'UseGeoCoords', true);
                    geoshow(land, 'FaceColor', [0.8 0.8 0.8]);hold on,alpha(0.7);
                    hold on,scatterm(SDtpFigS7(:,1)+0.5,SDtpFigS7(:,2)-0.5,5,SDtpFigS7(:,3));colormap(clr);clim([4 8]);%

tpslab=
save -ascii  SourceData_FigS7.txt SDtpFigS7;





%split and stats by regions
id=unique(cmat(:,3));nid=length(id);h = waitbar(0,'Please wait...');
lab={'EasternEurope';'NorthAmerica';'WesternEurope';'CentralAmerica';'SouthAmerica';'Australia';'Asia'};
for i=1:nid;waitbar(i/nid)
    fid=find(cmat(:,3)==id(i));smat=cmat(fid,4:end);
    out=quantile(smat,[0 0.025 0.25 0.5 0.75 0.975 1],1)
    out=[tps;out];eval(['save -ascii TempAnomalies_' cell2mat(lab(i)) '(stats).txt out']);
    clear out fid smat;    
end;close(h);

%%=========================================================================
%% STEP#3: FIGURES DISPLAY
%%=========================================================================
clear all, close all, clc,  
Weurope=load('TempAnomalies_WesternEurope(stats).txt');Eeurope=load('TempAnomalies_EasternEurope(stats).txt');
Namerica=load('TempAnomalies_NorthAmerica(stats).txt');Samerica=load('TempAnomalies_SouthAmerica(stats).txt');
Camerica=load('TempAnomalies_CentralAmerica(stats).txt');Australia=load('TempAnomalies_Australia(stats).txt');
Asia=load('TempAnomalies_Asia(stats).txt');

Weurope=Weurope-Weurope(:,1);Eeurope=Eeurope-Eeurope(:,1);Namerica=Namerica-Namerica(:,1);
Samerica=Samerica-Samerica(:,1);Camerica=Camerica-Camerica(:,1);Australia=Australia-Australia(:,1);
Asia=Asia-Asia(:,1);

cltrend=brewermap(8,'Set2');bsv=-100;clr1=brewermap(13,'BuPu');
figure,area(-Weurope(1,:),Weurope(7,:),'Facecolor',clr1(3,:),'Edgecolor','none','basevalue',bsv),
       hold on,area(-Weurope(1,:),Weurope(6,:),'Facecolor',clr1(5,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Weurope(1,:),Weurope(4,:),'Facecolor',clr1(3,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Weurope(1,:),Weurope(3,:),'Facecolor','w','Edgecolor','none','basevalue',bsv),
       hold on,plot(-Weurope(1,:),Weurope(5,:),'color',cltrend(3,:),'Linewidth',2),ylim([-25 5]),axis square;
       ylab('Temperature anomaly - Western Europe',13),xlab('Time (year BP)',13);
       tr=title('Late Pleistocene Temperature - Western Europe');set(tr,'Fontsize',14)

clr2=brewermap(13,'BuGn');
figure,area(-Eeurope(1,:),Eeurope(7,:),'Facecolor',clr2(5,:),'Edgecolor','none','basevalue',bsv),
       hold on,area(-Eeurope(1,:),Eeurope(6,:),'Facecolor',clr2(6,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Eeurope(1,:),Eeurope(4,:),'Facecolor',clr2(5,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Eeurope(1,:),Eeurope(3,:),'Facecolor','w','Edgecolor','none','basevalue',bsv),
       hold on,plot(-Eeurope(1,:),Eeurope(5,:),'color',cltrend(1,:),'Linewidth',2),ylim([-20 5]),axis square;
       ylab('Temperature anomaly - Eastern Europe',13),xlab('Time (year BP)',13);
       tr=title('Late Pleistocene Temperature - Eastern Europe');set(tr,'Fontsize',14)
       
clr3=brewermap(13,'RdGy');[~,c]=size(Namerica);for i=1:c;Namerica(2:end,i)=sort(Namerica(2:end,i),'ascend');end;
figure,area(-Namerica(1,:),Namerica(7,:),'Facecolor',clr3(6,:),'Edgecolor','none','basevalue',bsv),
       hold on,area(-Namerica(1,:),Namerica(6,:),'Facecolor',clr3(5,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Namerica(1,:),Namerica(4,:),'Facecolor',clr3(6,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Namerica(1,:),Namerica(3,:),'Facecolor','w','Edgecolor','none','basevalue',bsv),
       hold on,plot(-Namerica(1,:),Namerica(5,:),'color',cltrend(2,:),'Linewidth',2),ylim([-25 5]),axis square;
       ylab('Temperature anomaly - North America',13),xlab('Time (year BP)',13);
       tr=title('Late Pleistocene Temperature - North America');set(tr,'Fontsize',14)       

clr4=brewermap(13,'YlGn');[~,c]=size(Samerica);for i=1:c;Samerica(2:end,i)=sort(Samerica(2:end,i),'ascend');end;
figure,area(-Samerica(1,:),Samerica2(7,:),'Facecolor',clr4(4,:),'Edgecolor','none','basevalue',bsv),
       hold on,area(-Samerica(1,:),Samerica(6,:),'Facecolor',clr4(5,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Samerica(1,:),Samerica(4,:),'Facecolor',clr4(4,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Samerica(1,:),Samerica(3,:),'Facecolor','w','Edgecolor','none','basevalue',bsv),
       hold on,plot(-Samerica(1,:),Samerica(5,:),'color',cltrend(5,:),'Linewidth',2),ylim([-7 1]),axis square;
       ylab('Temperature anomaly - South America',13),xlab('Time (year BP)',13);
       tr=title('Late Pleistocene Temperature - South America');set(tr,'Fontsize',14)       

clr5=brewermap(13,'RdPu');[~,c]=size(Camerica);for i=1:c;Camerica(2:end,i)=sort(Camerica(2:end,i),'ascend');end;
figure,area(-Camerica(1,:),Camerica(7,:),'Facecolor',clr5(2,:),'Edgecolor','none','basevalue',bsv),
       hold on,area(-Camerica(1,:),Camerica(6,:),'Facecolor',clr5(4,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Camerica(1,:),Camerica(4,:),'Facecolor',clr5(2,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Camerica(1,:),Camerica(3,:),'Facecolor','w','Edgecolor','none','basevalue',bsv),
       hold on,plot(-Camerica(1,:),Camerica(5,:),'color',cltrend(4,:),'Linewidth',2),ylim([-5 1]),axis square;
       ylab('Temperature anomaly - Central America',13),xlab('Time (year BP)',13);
       tr=title('Late Pleistocene Temperature - Central America');set(tr,'Fontsize',14)       
       
clr6=brewermap(13,'YlOrBr');[~,c]=size(Camerica);for i=1:c;Australia(2:end,i)=sort(Australia(2:end,i),'ascend');end;
figure,area(-Australia(1,:),Australia(7,:),'Facecolor',clr6(3,:),'Edgecolor','none','basevalue',bsv),
       hold on,area(-Australia(1,:),Australia(6,:),'Facecolor',clr6(4,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Australia(1,:),Australia(4,:),'Facecolor',clr6(3,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Australia(1,:),Australia(3,:),'Facecolor','w','Edgecolor','none','basevalue',bsv),
       hold on,plot(-Australia(1,:),Australia(5,:),'color',cltrend(6,:),'Linewidth',2),ylim([-4 1]),axis square;
       ylab('Temperature anomaly - Australia',13),xlab('Time (year BP)',13);
       tr=title('Late Pleistocene Temperature - Australia');set(tr,'Fontsize',14)       
                           
clr7=brewermap(13,'BrBG');[~,c]=size(Asia);for i=1:c;Asia(2:end,i)=sort(Asia(2:end,i),'ascend');end;
figure,area(-Asia(1,:),Asia(7,:),'Facecolor',clr7(6,:),'Edgecolor','none','basevalue',bsv),
       hold on,area(-Asia(1,:),Asia(6,:),'Facecolor',clr7(5,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Asia(1,:),Asia(4,:),'Facecolor',clr7(6,:),'Edgecolor','none','basevalue',bsv)
       hold on,area(-Asia(1,:),Asia(3,:),'Facecolor','w','Edgecolor','none','basevalue',bsv),
       hold on,plot(-Asia(1,:),Asia(5,:),'color',clr7(3,:),'Linewidth',2),ylim([-9 1]),axis square;
       ylab('Temperature anomaly - Asia',13),xlab('Time (year BP)',13);
       tr=title('Late Pleistocene Temperature - Asia');set(tr,'Fontsize',14)     
       
       