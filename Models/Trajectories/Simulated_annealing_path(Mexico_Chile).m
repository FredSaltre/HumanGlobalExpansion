program Simulated_annealing_path(Mexico_Chile)
%on va calculer la vitesse et la direction de la colonoisation humaine
%le gradient temporel est la difference entre les valeur d'apparition
%calculees selon GRIWM et une dates la plus lointaoine qui est la premmiere
%date observeees dans les donnees

%%=========================================================================
%% VOIE #1
%%=========================================================================
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs=80;cs=299;re=143;ce=310;%startpt=max(max(mat.time));[rs,cs]=find(mat.time==startpt);[re,ce]=find(mat.time==endpt);endpt=min(min(mat.time));
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
[path, mat, ~]=shortpath(re,ce,rs,cs,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  

figure('Color','w'),worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs,cs),mat.lon(rs,cs),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re,ce),mat.lon(re,ce),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                   
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=2;MHent=[];TMHent=[];stoptime=20000;cpt=1;%0.12 est pas mal
figure(3),plot(cpt,Entropie,'ok','Markersize',4)%,hold on,plot(cpt,stop,'*b','Markersize',4);
figure(4),plot(cpt,TMH,'*b','Markersize',4);

while (stop < stoptime); %% on choisit un point au hasard et on le deplace dans son voisinage au hasard JUSQU'A AVOIR UN CERTAIN NOMBRE DE RUN CONSECUTIF SANS AMELIORATION
    megaflag=1;
    while (megaflag==1 && stop < stoptime) %%boucle fantoche qui n'a pour but de d'en sortir si on mne trouve pas une trajectoire satisfaisante => tombe dans la mer, etc...
        stop=stop+1
        pk=randi([path(2,1),path(end-1,1)],1);id=find(path(:,1)==pk);oldpath=path(id,:);
        cpt=cpt+1;
        %% on prend le voisin de ce point
        ngh=neighbourg(oldpath,mat);
        ngh(ngh(:,end)==1,:)=[];tfngh=isnan(ngh(:,5));ngh(tfngh==1,:)=[];
        if isempty(ngh==1);nghflag=1;break;end;
        pk2=randi(length(ngh(:,1)),1);newpt=ngh(pk2,:);clear ngh tfngh;
    
        path2=path;path2(id,2:end)=[newpt,NaN];path2(id,7)=1;%on change le parcours d'un au hasard 
   
        nghnew=neighbourg(path2(id,:),mat);nghnew(nghnew(:,1)==oldpath(2) & nghnew(:,2)==oldpath(3),:)=[]; %on regarde combien il y a de voisin au nouveau pts (on ne compte pas le point de depart)
        b4af=[path2(id-1,2:end);path2(id+1,2:end)];b4af=[b4af,zeros(2,1)];for i=1:2;idn=find(nghnew(:,1)==b4af(i,1) & nghnew(:,2)==b4af(i,2));if ~isempty(idn>0);b4af(i,end)=1;end;clear idn;end;
        clear nghnew;
    
        if (b4af(1,end)*b4af(2,end)==0);%si on perd la connexion avec les proche voisin => il faut rajouter des points pour retablir la connexion
            if and(b4af(1,end)==0,b4af(2,end)==0);%si on perd juste 1 voisin
                [bit1,~,spflag]=shortpath(newpt(1),newpt(2),b4af(1,1),b4af(1,2),mat);if (spflag==1);break;end;
                [bit2,~]=shortpath(b4af(2,1),b4af(2,2),newpt(1),newpt(2),mat);if (spflag==1);break;end;

                newbit=[bit1(1:end-1,:);bit2];
                b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));newpath=[path(1:b4-1,:);newbit];  
                aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));newpath=[newpath;path(aft+1:end,:)];  
                clear bit1 bit2;   
            else %si on perd les deux voisin a la fois           
                idn2=find(b4af(:,end)==0);[path3,~,spflag]=shortpath(b4af(idn2,1),b4af(idn2,2),newpt(1),newpt(2),mat);clear idn2;if (spflag==1);break;end;
                switch b4af(2,end)
                    case 0 %si c'est le voisin d'apres qui est perdu
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4,:);path3];aft=find(path(:,2)==newpath(end,2) & path(:,3)==newpath(end,3));
                        newpath=[newpath;path(aft+1:end,:)];
                    otherwise %si c'est le voisin d'avant qui est perdu
                        path3=flipud(path3);
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4-1,:);path3];aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));
                        newpath=[newpath;path(aft:end,:)];
                end;clear path3;
            end;clear b4 aft b4af;
        else %si on ne perd pas la connexion avec les proches voisins
            newpath=path2; 
        end;
        %%on recalcule la nouvelle entropie 
        p=proba3(newpath);Entropie2=sum(p(:,end));clear p path2; 

        %%% Comparaison des entropie
        MH=Entropie2/Entropie;
        if (MH<1)%cas ou l'on ameliore
            path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
            figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
            figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
        elseif (MH==1)% cas ou les entropie locales sont les meme => on rejete directement
            mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end; 
        else% cas ou l'on ameliore pas forcement => on met une condition sur l'acceptation
            pMH=exp(-(Entropie2-Entropie)/(TMH));pAcc=rand(1);
            if pAcc>(1-pMH);%on accepte en fonction d'une probabilite dependant des entropie locale + temperature
                path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
                figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
                figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
            else % on rejete en fonction d'une probabilite dependant des entropie locale + temperature
                mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;
            end;
            %figure(5),hold on,plot(cpt,pMH,'*b','Markersize',4);
            clear pAcc pMH;
        end;
        TMH=TMH-(MHrate/cpt);if TMH<=0;TMH=0.000000001;end;TMHent=[TMHent;TMH];
        clear MH newpath pk pk2 oldpath id Entropie2;
    end;
    stop=stop+1,
end;
pfinal=proba3(path);pfinal=[0;pfinal];


figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs,cs),mat.lon(rs,cs),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re,ce),mat.lon(re,ce),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Mexico-Chile_(path#1).txt path;save -ascii Trajectory_Mexico-Chile_(enthropy#1).txt MHent;
save -ascii Trajectory_Mexico-Chile_(Temp#1).txt TMHent;save -ascii Trajectory_Mexico-Chile_(Finalproba#1).txt pfinal;
   

%%=========================================================================
%% VOIE #2
%%=========================================================================
%% Mexico -> Brazil -> Chile 
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=80;cs_1=299;re_1=116;ce_1=331;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_2=116;cs_2=331;re_2=143;ce_2=310;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
path=[path_1;path_2];path(1:length(path(:,1)))=(1:length(path(:,1)))';mat.flag=mat_1.flag + mat_2.flag;[row,col]=find(mat.flag>1);mat.flag(row(:),col(:))=1;
id=find(path(:,2)<22 & path(:,3)<185);path(id,:)=[];mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;clear out nout;%on corrige la trajectoire pour les valeur abherentes
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);

                    
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=2.6;MHent=[];TMHent=[];stoptime=20000;cpt=1;%0.12 est pas mal
figure(3),plot(cpt,Entropie,'ok','Markersize',4)%,hold on,plot(cpt,stop,'*b','Markersize',4);
figure(4),plot(cpt,TMH,'*b','Markersize',4);

while (stop < stoptime); %% on choisit un point au hasard et on le deplace dans son voisinage au hasard JUSQU'A AVOIR UN CERTAIN NOMBRE DE RUN CONSECUTIF SANS AMELIORATION
    megaflag=1;
    while (megaflag==1 && stop < stoptime) %%boucle fantoche qui n'a pour but de d'en sortir si on mne trouve pas une trajectoire satisfaisante => tombe dans la mer, etc...
        stop=stop+1
        pk=randi([path(2,1),path(end-1,1)],1);id=find(path(:,1)==pk);oldpath=path(id,:);
        cpt=cpt+1;
        %% on prend le voisin de ce point
        ngh=neighbourg(oldpath,mat);
        ngh(ngh(:,end)==1,:)=[];tfngh=isnan(ngh(:,5));ngh(tfngh==1,:)=[];
        if isempty(ngh==1);nghflag=1;break;end;
        pk2=randi(length(ngh(:,1)),1);newpt=ngh(pk2,:);clear ngh tfngh;
    
        path2=path;path2(id,2:end)=[newpt,NaN];path2(id,7)=1;%on change le parcours d'un au hasard 
   
        nghnew=neighbourg(path2(id,:),mat);nghnew(nghnew(:,1)==oldpath(2) & nghnew(:,2)==oldpath(3),:)=[]; %on regarde combien il y a de voisin au nouveau pts (on ne compte pas le point de depart)
        b4af=[path2(id-1,2:end);path2(id+1,2:end)];b4af=[b4af,zeros(2,1)];for i=1:2;idn=find(nghnew(:,1)==b4af(i,1) & nghnew(:,2)==b4af(i,2));if ~isempty(idn>0);b4af(i,end)=1;end;clear idn;end;
        clear nghnew;
    
        if (b4af(1,end)*b4af(2,end)==0);%si on perd la connexion avec les proche voisin => il faut rajouter des points pour retablir la connexion
            if and(b4af(1,end)==0,b4af(2,end)==0);%si on perd juste 1 voisin
                [bit1,~,spflag]=shortpath(newpt(1),newpt(2),b4af(1,1),b4af(1,2),mat);if (spflag==1);break;end;
                [bit2,~]=shortpath(b4af(2,1),b4af(2,2),newpt(1),newpt(2),mat);if (spflag==1);break;end;

                newbit=[bit1(1:end-1,:);bit2];
                b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));newpath=[path(1:b4-1,:);newbit];  
                aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));newpath=[newpath;path(aft+1:end,:)];  
                clear bit1 bit2;   
            else %si on perd les deux voisin a la fois           
                idn2=find(b4af(:,end)==0);[path3,~,spflag]=shortpath(b4af(idn2,1),b4af(idn2,2),newpt(1),newpt(2),mat);clear idn2;if (spflag==1);break;end;
                switch b4af(2,end)
                    case 0 %si c'est le voisin d'apres qui est perdu
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4,:);path3];aft=find(path(:,2)==newpath(end,2) & path(:,3)==newpath(end,3));
                        newpath=[newpath;path(aft+1:end,:)];
                    otherwise %si c'est le voisin d'avant qui est perdu
                        path3=flipud(path3);
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4-1,:);path3];aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));
                        newpath=[newpath;path(aft:end,:)];
                end;clear path3;
            end;clear b4 aft b4af;
        else %si on ne perd pas la connexion avec les proches voisins
            newpath=path2; 
        end;
        %%on recalcule la nouvelle entropie 
        p=proba3(newpath);Entropie2=sum(p(:,end));clear p path2; 

        %%% Comparaison des entropie
        MH=Entropie2/Entropie;
        if (MH<1)%cas ou l'on ameliore
            path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
            figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
            figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
        elseif (MH==1)% cas ou les entropie locales sont les meme => on rejete directement
            mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end; 
        else% cas ou l'on ameliore pas forcement => on met une condition sur l'acceptation
            pMH=exp(-(Entropie2-Entropie)/(TMH));pAcc=rand(1);
            if pAcc>(1-pMH);%on accepte en fonction d'une probabilite dependant des entropie locale + temperature
                path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
                figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
                figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
            else % on rejete en fonction d'une probabilite dependant des entropie locale + temperature
                mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;
            end;
            %figure(5),hold on,plot(cpt,pMH,'*b','Markersize',4);
            clear pAcc pMH;
        end;
        TMH=TMH-(MHrate/cpt);if TMH<=0;TMH=0.000000001;end;TMHent=[TMHent;TMH];
        clear MH newpath pk pk2 oldpath id Entropie2;
    end;
    stop=stop+1,
end;
pfinal=proba3(path);pfinal=[0;pfinal];


figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Mexico-Chile_(path#2).txt path;save -ascii Trajectory_Mexico-Chile_(enthropy#2).txt MHent;
save -ascii Trajectory_Mexico-Chile_(Temp#2).txt TMHent;save -ascii Trajectory_Mexico-Chile_(Finalproba#2).txt pfinal;

%%=========================================================================
%% VOIE #3
%%=========================================================================
%% Mexico -> Brazil -> Chile 
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=80;cs_1=299;re_1=90;ce_1=329;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_2=90;cs_2=329;re_2=116;ce_2=331;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_3=116;cs_3=331;re_3=143;ce_3=310;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  

path=[path_1;path_2;path_3];path(1:length(path(:,1)))=(1:length(path(:,1)))';mat.flag=mat_1.flag + mat_2.flag + mat_3.flag;[row,col]=find(mat.flag>1);mat.flag(row(:),col(:))=1;
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);

                    
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=3.1;MHent=[];TMHent=[];stoptime=20000;cpt=1;%0.12 est pas mal
figure(3),plot(cpt,Entropie,'ok','Markersize',4)%,hold on,plot(cpt,stop,'*b','Markersize',4);
figure(4),plot(cpt,TMH,'*b','Markersize',4);

while (stop < stoptime); %% on choisit un point au hasard et on le deplace dans son voisinage au hasard JUSQU'A AVOIR UN CERTAIN NOMBRE DE RUN CONSECUTIF SANS AMELIORATION
    megaflag=1;
    while (megaflag==1 && stop < stoptime) %%boucle fantoche qui n'a pour but de d'en sortir si on mne trouve pas une trajectoire satisfaisante => tombe dans la mer, etc...
        stop=stop+1
        pk=randi([path(2,1),path(end-1,1)],1);id=find(path(:,1)==pk);oldpath=path(id,:);
        cpt=cpt+1;
        %% on prend le voisin de ce point
        ngh=neighbourg(oldpath,mat);
        ngh(ngh(:,end)==1,:)=[];tfngh=isnan(ngh(:,5));ngh(tfngh==1,:)=[];
        if isempty(ngh==1);nghflag=1;break;end;
        pk2=randi(length(ngh(:,1)),1);newpt=ngh(pk2,:);clear ngh tfngh;
    
        path2=path;path2(id,2:end)=[newpt,NaN];path2(id,7)=1;%on change le parcours d'un au hasard 
   
        nghnew=neighbourg(path2(id,:),mat);nghnew(nghnew(:,1)==oldpath(2) & nghnew(:,2)==oldpath(3),:)=[]; %on regarde combien il y a de voisin au nouveau pts (on ne compte pas le point de depart)
        b4af=[path2(id-1,2:end);path2(id+1,2:end)];b4af=[b4af,zeros(2,1)];for i=1:2;idn=find(nghnew(:,1)==b4af(i,1) & nghnew(:,2)==b4af(i,2));if ~isempty(idn>0);b4af(i,end)=1;end;clear idn;end;
        clear nghnew;
    
        if (b4af(1,end)*b4af(2,end)==0);%si on perd la connexion avec les proche voisin => il faut rajouter des points pour retablir la connexion
            if and(b4af(1,end)==0,b4af(2,end)==0);%si on perd juste 1 voisin
                [bit1,~,spflag]=shortpath(newpt(1),newpt(2),b4af(1,1),b4af(1,2),mat);if (spflag==1);break;end;
                [bit2,~]=shortpath(b4af(2,1),b4af(2,2),newpt(1),newpt(2),mat);if (spflag==1);break;end;

                newbit=[bit1(1:end-1,:);bit2];
                b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));newpath=[path(1:b4-1,:);newbit];  
                aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));newpath=[newpath;path(aft+1:end,:)];  
                clear bit1 bit2;   
            else %si on perd les deux voisin a la fois           
                idn2=find(b4af(:,end)==0);[path3,~,spflag]=shortpath(b4af(idn2,1),b4af(idn2,2),newpt(1),newpt(2),mat);clear idn2;if (spflag==1);break;end;
                switch b4af(2,end)
                    case 0 %si c'est le voisin d'apres qui est perdu
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4,:);path3];aft=find(path(:,2)==newpath(end,2) & path(:,3)==newpath(end,3));
                        newpath=[newpath;path(aft+1:end,:)];
                    otherwise %si c'est le voisin d'avant qui est perdu
                        path3=flipud(path3);
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4-1,:);path3];aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));
                        newpath=[newpath;path(aft:end,:)];
                end;clear path3;
            end;clear b4 aft b4af;
        else %si on ne perd pas la connexion avec les proches voisins
            newpath=path2; 
        end;
        %%on recalcule la nouvelle entropie 
        p=proba3(newpath);Entropie2=sum(p(:,end));clear p path2; 

        %%% Comparaison des entropie
        MH=Entropie2/Entropie;
        if (MH<1)%cas ou l'on ameliore
            path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
            figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
            figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
        elseif (MH==1)% cas ou les entropie locales sont les meme => on rejete directement
            mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end; 
        else% cas ou l'on ameliore pas forcement => on met une condition sur l'acceptation
            pMH=exp(-(Entropie2-Entropie)/(TMH));pAcc=rand(1);
            if pAcc>(1-pMH);%on accepte en fonction d'une probabilite dependant des entropie locale + temperature
                path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
                figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
                figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
            else % on rejete en fonction d'une probabilite dependant des entropie locale + temperature
                mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;
            end;
            %figure(5),hold on,plot(cpt,pMH,'*b','Markersize',4);
            clear pAcc pMH;
        end;
        TMH=TMH-(MHrate/cpt);if TMH<=0;TMH=0.000000001;end;TMHent=[TMHent;TMH];
        clear MH newpath pk pk2 oldpath id Entropie2;
    end;
    stop=stop+1,
end;
pfinal=proba3(path);pfinal=[0;pfinal];


figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Mexico-Chile_(path#3).txt path;save -ascii Trajectory_Mexico-Chile_(enthropy#3).txt MHent;
save -ascii Trajectory_Mexico-Chile_(Temp#3).txt TMHent;save -ascii Trajectory_Mexico-Chile_(Finalproba#3).txt pfinal;


%%=========================================================================
%% VOIE #4
%%=========================================================================
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=80;cs_1=299;re_1=90;ce_1=321;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_2=90;cs_2=321;re_2=112;ce_2=324;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_3=112;cs_3=324;re_3=143;ce_3=310;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
path=[path_1;path_2;path_3];path(1:length(path(:,1)))=(1:length(path(:,1)))';mat.flag=mat_1.flag + mat_2.flag + mat_3.flag;[row,col]=find(mat.flag>1);mat.flag(row(:),col(:))=1;
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_3,cs_3),mat.lon(rs_3,cs_3),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
    
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=3.8;MHent=[];TMHent=[];stoptime=20000;cpt=1;%0.12 est pas mal
figure(3),plot(cpt,Entropie,'ok','Markersize',4)%,hold on,plot(cpt,stop,'*b','Markersize',4);
figure(4),plot(cpt,TMH,'*b','Markersize',4);

while (stop < stoptime); %% on choisit un point au hasard et on le deplace dans son voisinage au hasard JUSQU'A AVOIR UN CERTAIN NOMBRE DE RUN CONSECUTIF SANS AMELIORATION
    megaflag=1;
    while (megaflag==1 && stop < stoptime) %%boucle fantoche qui n'a pour but de d'en sortir si on mne trouve pas une trajectoire satisfaisante => tombe dans la mer, etc...
        stop=stop+1
        pk=randi([path(2,1),path(end-1,1)],1);id=find(path(:,1)==pk);oldpath=path(id,:);
        cpt=cpt+1;
        %% on prend le voisin de ce point
        ngh=neighbourg(oldpath,mat);
        ngh(ngh(:,end)==1,:)=[];tfngh=isnan(ngh(:,5));ngh(tfngh==1,:)=[];
        if isempty(ngh==1);nghflag=1;break;end;
        pk2=randi(length(ngh(:,1)),1);newpt=ngh(pk2,:);clear ngh tfngh;
    
        path2=path;path2(id,2:end)=[newpt,NaN];path2(id,7)=1;%on change le parcours d'un au hasard 
   
        nghnew=neighbourg(path2(id,:),mat);nghnew(nghnew(:,1)==oldpath(2) & nghnew(:,2)==oldpath(3),:)=[]; %on regarde combien il y a de voisin au nouveau pts (on ne compte pas le point de depart)
        b4af=[path2(id-1,2:end);path2(id+1,2:end)];b4af=[b4af,zeros(2,1)];for i=1:2;idn=find(nghnew(:,1)==b4af(i,1) & nghnew(:,2)==b4af(i,2));if ~isempty(idn>0);b4af(i,end)=1;end;clear idn;end;
        clear nghnew;
    
        if (b4af(1,end)*b4af(2,end)==0);%si on perd la connexion avec les proche voisin => il faut rajouter des points pour retablir la connexion
            if and(b4af(1,end)==0,b4af(2,end)==0);%si on perd juste 1 voisin
                [bit1,~,spflag]=shortpath(newpt(1),newpt(2),b4af(1,1),b4af(1,2),mat);if (spflag==1);break;end;
                [bit2,~]=shortpath(b4af(2,1),b4af(2,2),newpt(1),newpt(2),mat);if (spflag==1);break;end;

                newbit=[bit1(1:end-1,:);bit2];
                b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));newpath=[path(1:b4-1,:);newbit];  
                aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));newpath=[newpath;path(aft+1:end,:)];  
                clear bit1 bit2;   
            else %si on perd les deux voisin a la fois           
                idn2=find(b4af(:,end)==0);[path3,~,spflag]=shortpath(b4af(idn2,1),b4af(idn2,2),newpt(1),newpt(2),mat);clear idn2;if (spflag==1);break;end;
                switch b4af(2,end)
                    case 0 %si c'est le voisin d'apres qui est perdu
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4,:);path3];aft=find(path(:,2)==newpath(end,2) & path(:,3)==newpath(end,3));
                        newpath=[newpath;path(aft+1:end,:)];
                    otherwise %si c'est le voisin d'avant qui est perdu
                        path3=flipud(path3);
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4-1,:);path3];aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));
                        newpath=[newpath;path(aft:end,:)];
                end;clear path3;
            end;clear b4 aft b4af;
        else %si on ne perd pas la connexion avec les proches voisins
            newpath=path2; 
        end;
        %%on recalcule la nouvelle entropie 
        p=proba3(newpath);Entropie2=sum(p(:,end));clear p path2; 

        %%% Comparaison des entropie
        MH=Entropie2/Entropie;
        if (MH<1)%cas ou l'on ameliore
            path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
            figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
            figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
        elseif (MH==1)% cas ou les entropie locales sont les meme => on rejete directement
            mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end; 
        else% cas ou l'on ameliore pas forcement => on met une condition sur l'acceptation
            pMH=exp(-(Entropie2-Entropie)/(TMH));pAcc=rand(1);
            if pAcc>(1-pMH);%on accepte en fonction d'une probabilite dependant des entropie locale + temperature
                path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
                figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
                figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
            else % on rejete en fonction d'une probabilite dependant des entropie locale + temperature
                mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;
            end;
            %figure(5),hold on,plot(cpt,pMH,'*b','Markersize',4);
            clear pAcc pMH;
        end;
        TMH=TMH-(MHrate/cpt);if TMH<=0;TMH=0.000000001;end;TMHent=[TMHent;TMH];
        clear MH newpath pk pk2 oldpath id Entropie2;
    end;
    stop=stop+1,
end;
pfinal=proba3(path);pfinal=[0;pfinal];


figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_3,cs_3),mat.lon(rs_3,cs_3),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Mexico-Chile_(path#4).txt path;save -ascii Trajectory_Mexico-Chile_(enthropy#4).txt MHent;
save -ascii Trajectory_Mexico-Chile_(Temp#4).txt TMHent;save -ascii Trajectory_Mexico-Chile_(Finalproba#4).txt pfinal;
   
%%=========================================================================
%% VOIE #5
%%=========================================================================
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=80;cs_1=299;re_1=90;ce_1=321;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_2=90;cs_2=321;re_2=112;ce_2=324;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_3=112;cs_3=324;re_3=129;ce_3=324;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_4=129;cs_4=324;re_4=136;ce_4=317;[path_4, mat_4, ~]=shortpath(re_4,ce_4,rs_4,cs_4,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_5=136;cs_5=317;re_5=143;ce_5=310;[path_5, mat_5, ~]=shortpath(re_5,ce_5,rs_5,cs_5,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
path=[path_1;path_2;path_3;path_4;path_5];path(1:length(path(:,1)))=(1:length(path(:,1)))';mat.flag=mat_1.flag + mat_2.flag + mat_3.flag + mat_4.flag + mat_5.flag;[row,col]=find(mat.flag>1);mat.flag(row(:),col(:))=1;
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_5,cs_5),mat.lon(rs_5,cs_5),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
    
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=4.2;MHent=[];TMHent=[];stoptime=20000;cpt=1;%0.12 est pas mal
figure(3),plot(cpt,Entropie,'ok','Markersize',4)%,hold on,plot(cpt,stop,'*b','Markersize',4);
figure(4),plot(cpt,TMH,'*b','Markersize',4);

while (stop < stoptime); %% on choisit un point au hasard et on le deplace dans son voisinage au hasard JUSQU'A AVOIR UN CERTAIN NOMBRE DE RUN CONSECUTIF SANS AMELIORATION
    megaflag=1;
    while (megaflag==1 && stop < stoptime) %%boucle fantoche qui n'a pour but de d'en sortir si on mne trouve pas une trajectoire satisfaisante => tombe dans la mer, etc...
        stop=stop+1
        pk=randi([path(2,1),path(end-1,1)],1);id=find(path(:,1)==pk);oldpath=path(id,:);
        cpt=cpt+1;
        %% on prend le voisin de ce point
        ngh=neighbourg(oldpath,mat);
        ngh(ngh(:,end)==1,:)=[];tfngh=isnan(ngh(:,5));ngh(tfngh==1,:)=[];
        if isempty(ngh==1);nghflag=1;break;end;
        pk2=randi(length(ngh(:,1)),1);newpt=ngh(pk2,:);clear ngh tfngh;
    
        path2=path;path2(id,2:end)=[newpt,NaN];path2(id,7)=1;%on change le parcours d'un au hasard 
   
        nghnew=neighbourg(path2(id,:),mat);nghnew(nghnew(:,1)==oldpath(2) & nghnew(:,2)==oldpath(3),:)=[]; %on regarde combien il y a de voisin au nouveau pts (on ne compte pas le point de depart)
        b4af=[path2(id-1,2:end);path2(id+1,2:end)];b4af=[b4af,zeros(2,1)];for i=1:2;idn=find(nghnew(:,1)==b4af(i,1) & nghnew(:,2)==b4af(i,2));if ~isempty(idn>0);b4af(i,end)=1;end;clear idn;end;
        clear nghnew;
    
        if (b4af(1,end)*b4af(2,end)==0);%si on perd la connexion avec les proche voisin => il faut rajouter des points pour retablir la connexion
            if and(b4af(1,end)==0,b4af(2,end)==0);%si on perd juste 1 voisin
                [bit1,~,spflag]=shortpath(newpt(1),newpt(2),b4af(1,1),b4af(1,2),mat);if (spflag==1);break;end;
                [bit2,~]=shortpath(b4af(2,1),b4af(2,2),newpt(1),newpt(2),mat);if (spflag==1);break;end;

                newbit=[bit1(1:end-1,:);bit2];
                b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));newpath=[path(1:b4-1,:);newbit];  
                aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));newpath=[newpath;path(aft+1:end,:)];  
                clear bit1 bit2;   
            else %si on perd les deux voisin a la fois           
                idn2=find(b4af(:,end)==0);[path3,~,spflag]=shortpath(b4af(idn2,1),b4af(idn2,2),newpt(1),newpt(2),mat);clear idn2;if (spflag==1);break;end;
                switch b4af(2,end)
                    case 0 %si c'est le voisin d'apres qui est perdu
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4,:);path3];aft=find(path(:,2)==newpath(end,2) & path(:,3)==newpath(end,3));
                        newpath=[newpath;path(aft+1:end,:)];
                    otherwise %si c'est le voisin d'avant qui est perdu
                        path3=flipud(path3);
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4-1,:);path3];aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));
                        newpath=[newpath;path(aft:end,:)];
                end;clear path3;
            end;clear b4 aft b4af;
        else %si on ne perd pas la connexion avec les proches voisins
            newpath=path2; 
        end;
        %%on recalcule la nouvelle entropie 
        p=proba3(newpath);Entropie2=sum(p(:,end));clear p path2; 

        %%% Comparaison des entropie
        MH=Entropie2/Entropie;
        if (MH<1)%cas ou l'on ameliore
            path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
            figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
            figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
        elseif (MH==1)% cas ou les entropie locales sont les meme => on rejete directement
            mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end; 
        else% cas ou l'on ameliore pas forcement => on met une condition sur l'acceptation
            pMH=exp(-(Entropie2-Entropie)/(TMH));pAcc=rand(1);
            if pAcc>(1-pMH);%on accepte en fonction d'une probabilite dependant des entropie locale + temperature
                path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
                figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
                figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
            else % on rejete en fonction d'une probabilite dependant des entropie locale + temperature
                mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;
            end;
            %figure(5),hold on,plot(cpt,pMH,'*b','Markersize',4);
            clear pAcc pMH;
        end;
        TMH=TMH-(MHrate/cpt);if TMH<=0;TMH=0.000000001;end;TMHent=[TMHent;TMH];
        clear MH newpath pk pk2 oldpath id Entropie2;
    end;
    stop=stop+1,
end;
pfinal=proba3(path);pfinal=[0;pfinal];

figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_5,cs_5),mat.lon(rs_5,cs_5),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Mexico-Chile_(path#5).txt path;save -ascii Trajectory_Mexico-Chile_(enthropy#5).txt MHent;
save -ascii Trajectory_Mexico-Chile_(Temp#5).txt TMHent;save -ascii Trajectory_Mexico-Chile_(Finalproba#5).txt pfinal;
   
%%=========================================================================
%% VOIE #6
%%=========================================================================
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=80;cs_1=299;re_1=93;ce_1=301;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_2=93;cs_2=301;re_2=106;ce_2=305;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_3=106;cs_3=305;re_3=108;ce_3=310;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_4=108;cs_4=310;re_4=118;ce_4=314;[path_4, mat_4, ~]=shortpath(re_4,ce_4,rs_4,cs_4,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_5=118;cs_5=314;re_5=132;ce_5=309;[path_5, mat_5, ~]=shortpath(re_5,ce_5,rs_5,cs_5,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_6=132;cs_6=309;re_6=143;ce_6=310;[path_6, mat_6, ~]=shortpath(re_6,ce_6,rs_6,cs_6,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
path=[path_1;path_2;path_3;path_4;path_5;path_6];path(1:length(path(:,1)))=(1:length(path(:,1)))';mat.flag=mat_1.flag + mat_2.flag + mat_3.flag + mat_4.flag + mat_5.flag + mat_6.flag;[row,col]=find(mat.flag>1);mat.flag(row(:),col(:))=1;
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_6,cs_6),mat.lon(rs_6,cs_6),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
    
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=4.35;MHent=[];TMHent=[];stoptime=20000;cpt=1;%0.12 est pas mal
figure(3),plot(cpt,Entropie,'ok','Markersize',4)%,hold on,plot(cpt,stop,'*b','Markersize',4);
figure(4),plot(cpt,TMH,'*b','Markersize',4);

while (stop < stoptime); %% on choisit un point au hasard et on le deplace dans son voisinage au hasard JUSQU'A AVOIR UN CERTAIN NOMBRE DE RUN CONSECUTIF SANS AMELIORATION
    megaflag=1;
    while (megaflag==1 && stop < stoptime) %%boucle fantoche qui n'a pour but de d'en sortir si on mne trouve pas une trajectoire satisfaisante => tombe dans la mer, etc...
        stop=stop+1
        pk=randi([path(2,1),path(end-1,1)],1);id=find(path(:,1)==pk);oldpath=path(id,:);
        cpt=cpt+1;
        %% on prend le voisin de ce point
        ngh=neighbourg(oldpath,mat);
        ngh(ngh(:,end)==1,:)=[];tfngh=isnan(ngh(:,5));ngh(tfngh==1,:)=[];
        if isempty(ngh==1);nghflag=1;break;end;
        pk2=randi(length(ngh(:,1)),1);newpt=ngh(pk2,:);clear ngh tfngh;
    
        path2=path;path2(id,2:end)=[newpt,NaN];path2(id,7)=1;%on change le parcours d'un au hasard 
   
        nghnew=neighbourg(path2(id,:),mat);nghnew(nghnew(:,1)==oldpath(2) & nghnew(:,2)==oldpath(3),:)=[]; %on regarde combien il y a de voisin au nouveau pts (on ne compte pas le point de depart)
        b4af=[path2(id-1,2:end);path2(id+1,2:end)];b4af=[b4af,zeros(2,1)];for i=1:2;idn=find(nghnew(:,1)==b4af(i,1) & nghnew(:,2)==b4af(i,2));if ~isempty(idn>0);b4af(i,end)=1;end;clear idn;end;
        clear nghnew;
    
        if (b4af(1,end)*b4af(2,end)==0);%si on perd la connexion avec les proche voisin => il faut rajouter des points pour retablir la connexion
            if and(b4af(1,end)==0,b4af(2,end)==0);%si on perd juste 1 voisin
                [bit1,~,spflag]=shortpath(newpt(1),newpt(2),b4af(1,1),b4af(1,2),mat);if (spflag==1);break;end;
                [bit2,~]=shortpath(b4af(2,1),b4af(2,2),newpt(1),newpt(2),mat);if (spflag==1);break;end;

                newbit=[bit1(1:end-1,:);bit2];
                b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));newpath=[path(1:b4-1,:);newbit];  
                aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));newpath=[newpath;path(aft+1:end,:)];  
                clear bit1 bit2;   
            else %si on perd les deux voisin a la fois           
                idn2=find(b4af(:,end)==0);[path3,~,spflag]=shortpath(b4af(idn2,1),b4af(idn2,2),newpt(1),newpt(2),mat);clear idn2;if (spflag==1);break;end;
                switch b4af(2,end)
                    case 0 %si c'est le voisin d'apres qui est perdu
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4,:);path3];aft=find(path(:,2)==newpath(end,2) & path(:,3)==newpath(end,3));
                        newpath=[newpath;path(aft+1:end,:)];
                    otherwise %si c'est le voisin d'avant qui est perdu
                        path3=flipud(path3);
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4-1,:);path3];aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));
                        newpath=[newpath;path(aft:end,:)];
                end;clear path3;
            end;clear b4 aft b4af;
        else %si on ne perd pas la connexion avec les proches voisins
            newpath=path2; 
        end;
        %%on recalcule la nouvelle entropie 
        p=proba3(newpath);Entropie2=sum(p(:,end));clear p path2; 

        %%% Comparaison des entropie
        MH=Entropie2/Entropie;
        if (MH<1)%cas ou l'on ameliore
            path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
            figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
            figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
        elseif (MH==1)% cas ou les entropie locales sont les meme => on rejete directement
            mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end; 
        else% cas ou l'on ameliore pas forcement => on met une condition sur l'acceptation
            pMH=exp(-(Entropie2-Entropie)/(TMH));pAcc=rand(1);
            if pAcc>(1-pMH);%on accepte en fonction d'une probabilite dependant des entropie locale + temperature
                path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
                figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
                figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
            else % on rejete en fonction d'une probabilite dependant des entropie locale + temperature
                mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;
            end;
            %figure(5),hold on,plot(cpt,pMH,'*b','Markersize',4);
            clear pAcc pMH;
        end;
        TMH=TMH-(MHrate/cpt);if TMH<=0;TMH=0.000000001;end;TMHent=[TMHent;TMH];
        clear MH newpath pk pk2 oldpath id Entropie2;
    end;
    stop=stop+1,
end;
pfinal=proba3(path);pfinal=[0;pfinal];

figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_5,cs_5),mat.lon(rs_5,cs_5),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Mexico-Chile_(path#6).txt path;save -ascii Trajectory_Mexico-Chile_(enthropy#6).txt MHent;
save -ascii Trajectory_Mexico-Chile_(Temp#6).txt TMHent;save -ascii Trajectory_Mexico-Chile_(Finalproba#6).txt pfinal;

%%=========================================================================
%% VOIE #7
%%=========================================================================
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=80;cs_1=299;re_1=94;ce_1=340;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_2=94;cs_2=340;re_2=108;ce_2=342;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_3=108;cs_3=342;re_3=115;ce_3=334;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_4=115;cs_4=334;re_4=143;ce_4=310;[path_4, mat_4, ~]=shortpath(re_4,ce_4,rs_4,cs_4,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
path=[path_1;path_2;path_3;path_4];path(1:length(path(:,1)))=(1:length(path(:,1)))';mat.flag=mat_1.flag + mat_2.flag + mat_3.flag + mat_4.flag;[row,col]=find(mat.flag>1);mat.flag(row(:),col(:))=1;
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_4,cs_4),mat.lon(rs_4,cs_4),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
    
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=3.9;MHent=[];TMHent=[];stoptime=20000;cpt=1;%0.12 est pas mal
figure(3),plot(cpt,Entropie,'ok','Markersize',4)%,hold on,plot(cpt,stop,'*b','Markersize',4);
figure(4),plot(cpt,TMH,'*b','Markersize',4);

while (stop < stoptime); %% on choisit un point au hasard et on le deplace dans son voisinage au hasard JUSQU'A AVOIR UN CERTAIN NOMBRE DE RUN CONSECUTIF SANS AMELIORATION
    megaflag=1;
    while (megaflag==1 && stop < stoptime) %%boucle fantoche qui n'a pour but de d'en sortir si on mne trouve pas une trajectoire satisfaisante => tombe dans la mer, etc...
        stop=stop+1
        pk=randi([path(2,1),path(end-1,1)],1);id=find(path(:,1)==pk);oldpath=path(id,:);
        cpt=cpt+1;
        %% on prend le voisin de ce point
        ngh=neighbourg(oldpath,mat);
        ngh(ngh(:,end)==1,:)=[];tfngh=isnan(ngh(:,5));ngh(tfngh==1,:)=[];
        if isempty(ngh==1);nghflag=1;break;end;
        pk2=randi(length(ngh(:,1)),1);newpt=ngh(pk2,:);clear ngh tfngh;
    
        path2=path;path2(id,2:end)=[newpt,NaN];path2(id,7)=1;%on change le parcours d'un au hasard 
   
        nghnew=neighbourg(path2(id,:),mat);nghnew(nghnew(:,1)==oldpath(2) & nghnew(:,2)==oldpath(3),:)=[]; %on regarde combien il y a de voisin au nouveau pts (on ne compte pas le point de depart)
        b4af=[path2(id-1,2:end);path2(id+1,2:end)];b4af=[b4af,zeros(2,1)];for i=1:2;idn=find(nghnew(:,1)==b4af(i,1) & nghnew(:,2)==b4af(i,2));if ~isempty(idn>0);b4af(i,end)=1;end;clear idn;end;
        clear nghnew;
    
        if (b4af(1,end)*b4af(2,end)==0);%si on perd la connexion avec les proche voisin => il faut rajouter des points pour retablir la connexion
            if and(b4af(1,end)==0,b4af(2,end)==0);%si on perd juste 1 voisin
                [bit1,~,spflag]=shortpath(newpt(1),newpt(2),b4af(1,1),b4af(1,2),mat);if (spflag==1);break;end;
                [bit2,~]=shortpath(b4af(2,1),b4af(2,2),newpt(1),newpt(2),mat);if (spflag==1);break;end;

                newbit=[bit1(1:end-1,:);bit2];
                b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));newpath=[path(1:b4-1,:);newbit];  
                aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));newpath=[newpath;path(aft+1:end,:)];  
                clear bit1 bit2;   
            else %si on perd les deux voisin a la fois           
                idn2=find(b4af(:,end)==0);[path3,~,spflag]=shortpath(b4af(idn2,1),b4af(idn2,2),newpt(1),newpt(2),mat);clear idn2;if (spflag==1);break;end;
                switch b4af(2,end)
                    case 0 %si c'est le voisin d'apres qui est perdu
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4,:);path3];aft=find(path(:,2)==newpath(end,2) & path(:,3)==newpath(end,3));
                        newpath=[newpath;path(aft+1:end,:)];
                    otherwise %si c'est le voisin d'avant qui est perdu
                        path3=flipud(path3);
                        b4=find(path(:,2)==b4af(1,1) & path(:,3)==b4af(1,2));
                        newpath=[path(1:b4-1,:);path3];aft=find(path(:,2)==b4af(2,1) & path(:,3)==b4af(2,2));
                        newpath=[newpath;path(aft:end,:)];
                end;clear path3;
            end;clear b4 aft b4af;
        else %si on ne perd pas la connexion avec les proches voisins
            newpath=path2; 
        end;
        %%on recalcule la nouvelle entropie 
        p=proba3(newpath);Entropie2=sum(p(:,end));clear p path2; 

        %%% Comparaison des entropie
        MH=Entropie2/Entropie;
        if (MH<1)%cas ou l'on ameliore
            path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
            figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
            figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
        elseif (MH==1)% cas ou les entropie locales sont les meme => on rejete directement
            mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end; 
        else% cas ou l'on ameliore pas forcement => on met une condition sur l'acceptation
            pMH=exp(-(Entropie2-Entropie)/(TMH));pAcc=rand(1);
            if pAcc>(1-pMH);%on accepte en fonction d'une probabilite dependant des entropie locale + temperature
                path=newpath;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;Entropie=Entropie2;stop=0;MHent=[MHent;Entropie];
                figure(3),hold on,plot(cpt,Entropie,'ok','Markersize',4),%hold on,plot(cpt,stop,'*b','Markersize',4);
                figure(4),hold on,plot(cpt,TMH,'*b','Markersize',4);
            else % on rejete en fonction d'une probabilite dependant des entropie locale + temperature
                mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;
            end;
            %figure(5),hold on,plot(cpt,pMH,'*b','Markersize',4);
            clear pAcc pMH;
        end;
        TMH=TMH-(MHrate/cpt);if TMH<=0;TMH=0.000000001;end;TMHent=[TMHent;TMH];
        clear MH newpath pk pk2 oldpath id Entropie2;
    end;
    stop=stop+1,
end;
pfinal=proba3(path);pfinal=[0;pfinal];

figure('Color','w'),ax=worldmap({'North America','South America'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_4,cs_4),mat.lon(rs_4,cs_4),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Mexico-Chile_(path#7).txt path;save -ascii Trajectory_Mexico-Chile_(enthropy#7).txt MHent;
save -ascii Trajectory_Mexico-Chile_(Temp#7).txt TMHent;save -ascii Trajectory_Mexico-Chile_(Finalproba#7).txt pfinal;

%%=========================================================================
%% SHORT PATHS
%%=========================================================================
% compilation de toutes les trajectoires initiales avant minimisation
% short path #1 => the shortest to use
clear all, close all, clc,
mat.time=load('Humancolonisation_timing_v4.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs=80;cs=299;re=143;ce=310;;[path, mat, ~]=shortpath(re,ce,rs,cs,mat);
out=[];for i=90:109;for j=300:310;out=[out;i,j];end;end;nout=length(out(:,1));for i=1:nout;id=find(path(:,2)==out(i,1) & path(:,3)==out(i,2));if length(id)>=1; path(id,:)=[];end;end;mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;clear out nout;
p=proba3(path);Entropie=sum(p);
save -ascii Short_Trajectory_Mexico-Chile_(path#1).txt path;save -ascii Short_Trajectory_Mexico-Chile_(enthropy#1).txt Entropie;save -ascii Short_Trajectory_Mexico-Chile_(Finalproba#1).txt p;

% short path #2  
clear all, close all, clc,
mat.time=load('Humancolonisation_timing_v4.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=80;cs_1=299;re_1=116;ce_1=331;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);rs_2=116;cs_2=331;re_2=143;ce_2=310;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);
path=[path_1;path_2];id=find(path(:,2)<22 & path(:,3)<185);path(id,:)=[];mat.flag(:,:)=0;for i=1:length(path(:,1));mat.flag(path(i,2),path(i,3))=1;path(i,1)=i;end;clear out nout;%on corrige la trajectoire pour les valeur abherentes
p=proba2(path);Entropie=sum(p);
save -ascii Short_Trajectory_Mexico-Chile_(path#2).txt path;save -ascii Short_Trajectory_Mexico-Chile_(enthropy#2).txt Entropie;save -ascii Short_Trajectory_Mexico-Chile_(Finalproba#2).txt p;

% short path #3  
clear all, close all, clc,
mat.time=load('Humancolonisation_timing_v4.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=80;cs_1=299;re_1=90;ce_1=329;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);rs_2=90;cs_2=329;re_2=116;ce_2=331;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);
rs_3=116;cs_3=331;re_3=143;ce_3=310;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);path=[path_1;path_2;path_3];p=proba2(path);Entropie=sum(p);
save -ascii Short_Trajectory_Mexico-Chile_(path#3).txt path;save -ascii Short_Trajectory_Mexico-Chile_(enthropy#3).txt Entropie;save -ascii Short_Trajectory_Mexico-Chile_(Finalproba#3).txt p;

% short path #4  
clear all, close all, clc,
mat.time=load('Humancolonisation_timing_v4.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=80;cs_1=299;re_1=90;ce_1=321;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);rs_2=90;cs_2=321;re_2=112;ce_2=324;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);
rs_3=112;cs_3=324;re_3=143;ce_3=310;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);path=[path_1;path_2;path_3];p=proba2(path);Entropie=sum(p);
save -ascii Short_Trajectory_Mexico-Chile_(path#4).txt path;save -ascii Short_Trajectory_Mexico-Chile_(enthropy#4).txt Entropie;save -ascii Short_Trajectory_Mexico-Chile_(Finalproba#4).txt p;

% short path #5  
clear all, close all, clc,
mat.time=load('Humancolonisation_timing_v4.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=80;cs_1=299;re_1=90;ce_1=321;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);rs_2=90;cs_2=321;re_2=112;ce_2=324;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);
rs_3=112;cs_3=324;re_3=129;ce_3=324;;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);rs_4=129;cs_4=324;re_4=136;ce_4=317;[path_4, mat_4, ~]=shortpath(re_4,ce_4,rs_4,cs_4,mat);%
rs_5=136;cs_5=317;re_5=143;ce_5=310;[path_5, mat_5, ~]=shortpath(re_5,ce_5,rs_5,cs_5,mat);
path=[path_1;path_2;path_3;path_4;path_5];p=proba2(path);Entropie=sum(p);
save -ascii Short_Trajectory_Mexico-Chile_(path#5).txt path;save -ascii Short_Trajectory_Mexico-Chile_(enthropy#5).txt Entropie;save -ascii Short_Trajectory_Mexico-Chile_(Finalproba#5).txt p;

% short path #6  
clear all, close all, clc,
mat.time=load('Humancolonisation_timing_v4.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=80;cs_1=299;re_1=93;ce_1=301;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);rs_2=93;cs_2=301;re_2=106;ce_2=305;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);
rs_3=106;cs_3=305;re_3=108;ce_3=310;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);rs_4=108;cs_4=310;re_4=118;ce_4=314;[path_4, mat_4, ~]=shortpath(re_4,ce_4,rs_4,cs_4,mat);
rs_5=118;cs_5=314;re_5=132;ce_5=309;[path_5, mat_5, ~]=shortpath(re_5,ce_5,rs_5,cs_5,mat);rs_6=132;cs_6=309;re_6=143;ce_6=310;[path_6, mat_6, ~]=shortpath(re_6,ce_6,rs_6,cs_6,mat);
path=[path_1;path_2;path_3;path_4;path_5;path_6];p=proba2(path);Entropie=sum(p);
save -ascii Short_Trajectory_Mexico-Chile_(path#6).txt path;save -ascii Short_Trajectory_Mexico-Chile_(enthropy#6).txt Entropie;save -ascii Short_Trajectory_Mexico-Chile_(Finalproba#6).txt p;

