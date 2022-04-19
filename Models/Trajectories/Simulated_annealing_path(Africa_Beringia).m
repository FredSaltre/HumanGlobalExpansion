program Simulated_annealing_path(Africa_Beringia)
%on va calculer la vitesse et la direction de la colonoisation humaine
%le gradient temporel est la difference entre les valeur d'apparition
%calculees selon GRIWM et une dates la plus lointaoine qui est la premmiere
%date observeees dans les donnees

% pour l'optimisation on joue sur deux parametres: TMH et MHrate
% ## plus MHrate est grand on plus on va durcir l'optimisation rapidement et
%    donc minimiser tres vite mais en exp[lorant peu l'espace. Plus il est
%    petit et plus l'algorithma va faire des recherche et eviter les minimum
%    locaux. Toutefois s'il est trop petit on peut tres bien par ne plus
%    arriver a minimiser l'enthropie car on accepte trop de valeur qui ne
%    minimisent pas 
% ## THM est la temperature de depart, c'est le meme principe que MHrate,
%    il ne fait pas qu'elle soit trop grande sinon on mets un tem,ps fou a
%    minimiser

%%=========================================================================
%% VOIE #1
%%=========================================================================
%% shortest path
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs=60;cs=62;re=24;ce=178;%startpt=max(max(mat.time));[rs,cs]=find(mat.time==startpt);[re,ce]=find(mat.time==endpt);endpt=min(min(mat.time));
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
[path, mat, ~]=shortpath(re,ce,rs,cs,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
figure('Color','w'),ax=worldmap({'Europe','Asia'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs,cs),mat.lon(rs,cs),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re,ce),mat.lon(re,ce),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                  
                    
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=2.2;MHent=[];TMHent=[];stoptime=10000;cpt=1;%0.12 est pas mal
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

figure('Color','w'),ax=worldmap({'Europe','Asia'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs,cs),mat.lon(rs,cs),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re,ce),mat.lon(re,ce),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);

figure('Color','w'),ax=worldmap('World');setm(ax, 'Origin', [0 160 0]);gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs,cs),mat.lon(rs,cs),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re,ce),mat.lon(re,ce),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Africa-Beringia_(path#1).txt path;save -ascii Trajectory_Africa-Beringia_(enthropy#1).txt MHent;
save -ascii Trajectory_Africa-Beringia_(Temp#1).txt TMHent;save -ascii Trajectory_Africa-Beringia_(Finalproba#1).txt pfinal;

%%=========================================================================
%% VOIE #2
%%=========================================================================
%% africa -> germany -> Beringia 
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=60;cs_1=62;re_1=48;ce_1=87;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_2=48;cs_2=87;re_2=24;ce_2=178;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
path=[path_1;path_2];path(1:length(path(:,1)))=(1:length(path(:,1)))';mat.flag=mat_1.flag + mat_2.flag;[row,col]=find(mat.flag>1);mat.flag(row(:),col(:))=1;
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
figure('Color','w'),ax=worldmap({'Europe','Asia'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                  
                    
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=3;MHent=[];TMHent=[];stoptime=10000;cpt=1;%0.12 est pas mal
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

figure('Color','w'),ax=worldmap({'Europe','Asia'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);

figure('Color','w'),ax=worldmap('World');setm(ax, 'Origin', [0 160 0]);gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs,cs),mat.lon(rs,cs),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re,ce),mat.lon(re,ce),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Africa-Beringia_(path#2).txt path;save -ascii Trajectory_Africa-Beringia_(enthropy#2).txt MHent;
save -ascii Trajectory_Africa-Beringia_(Temp#2).txt TMHent;save -ascii Trajectory_Africa-Beringia_(Finalproba#2).txt pfinal;

%%=========================================================================
%% VOIE #3
%%=========================================================================
%% africa -> under mer caspienne -> east mer caspienne -> finlande 
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=60;cs_1=62;re_1=44;ce_1=127;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_2=44;cs_2=127;re_2=24;ce_2=178;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
path=[path_1;path_2];path(1:length(path(:,1)))=(1:length(path(:,1)))';mat.flag=mat_1.flag + mat_2.flag;[row,col]=find(mat.flag>1);mat.flag(row(:),col(:))=1;
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
figure('Color','w'),ax=worldmap({'Europe','Asia'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                  
                    
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=2.2;MHent=[];TMHent=[];stoptime=20000;cpt=1;%0.12 est pas mal
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

figure('Color','w'),ax=worldmap({'Europe','Asia'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);

figure('Color','w'),ax=worldmap('World');setm(ax, 'Origin', [0 160 0]);gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs,cs),mat.lon(rs,cs),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re,ce),mat.lon(re,ce),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Africa-Beringia_(path#3).txt path;save -ascii Trajectory_Africa-Beringia_(enthropy#3).txt MHent;
save -ascii Trajectory_Africa-Beringia_(Temp#3).txt TMHent;save -ascii Trajectory_Africa-Beringia_(Finalproba#3).txt pfinal;

%%=========================================================================
%% VOIE #4
%%=========================================================================
%% africa -> North Australia -> Beringia 
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=60;cs_1=62;re_1=76;ce_1=125;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_2=76;cs_2=125;re_2=24;ce_2=178;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
path=[path_1;path_2];path(1:length(path(:,1)))=(1:length(path(:,1)))';mat.flag=mat_1.flag + mat_2.flag;[row,col]=find(mat.flag>1);mat.flag(row(:),col(:))=1;
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
figure('Color','w'),ax=worldmap({'Europe','Asia'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                  
                    
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=3.3;MHent=[];TMHent=[];stoptime=20000;cpt=1;%0.12 est pas mal
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

figure('Color','w'),ax=worldmap({'Europe','Asia'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);

figure('Color','w'),ax=worldmap('World');setm(ax, 'Origin', [0 160 0]);gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs,cs),mat.lon(rs,cs),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re,ce),mat.lon(re,ce),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Africa-Beringia_(path#4).txt path;save -ascii Trajectory_Africa-Beringia_(enthropy#4).txt MHent;
save -ascii Trajectory_Africa-Beringia_(Temp#4).txt TMHent;save -ascii Trajectory_Africa-Beringia_(Finalproba#4).txt pfinal;

%%=========================================================================
%% VOIE #5
%%=========================================================================
%% africa -> Chine -> Beringia 
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=60;cs_1=62;re_1=58;ce_1=114;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_2=58;cs_2=114;re_2=46;ce_2=135;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_3=46;cs_3=135;re_3=24;ce_3=178;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
path=[path_1;path_2;path_3];path(1:length(path(:,1)))=(1:length(path(:,1)))';mat.flag=mat_1.flag+mat_2.flag+mat_3.flag;[row,col]=find(mat.flag>1);mat.flag(row(:),col(:))=1;
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
figure('Color','w'),ax=worldmap({'Europe','Asia'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                  
                    
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=2.1;MHent=[];TMHent=[];stoptime=20000;cpt=1;%0.12 est pas mal
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


figure('Color','w'),ax=worldmap({'Europe','Asia'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);

figure('Color','w'),ax=worldmap('World');setm(ax, 'Origin', [0 160 0]);gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs,cs),mat.lon(rs,cs),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re,ce),mat.lon(re,ce),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Africa-Beringia_(path#5).txt path;save -ascii Trajectory_Africa-Beringia_(enthropy#5).txt MHent;
save -ascii Trajectory_Africa-Beringia_(Temp#5).txt TMHent;save -ascii Trajectory_Africa-Beringia_(Finalproba#5).txt pfinal;

%%=========================================================================
%% VOIE #6
%%=========================================================================
%% africa -> North France -> Poland --> Beringia 
clear all, close all, clc,
mat.time=load('maskedHumantimingraster_v2.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);%st=50;mat.time(:,st:end)=NaN;mat.time(17:20,34:38)=NaN;mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=60;cs_1=62;re_1=44;ce_1=25;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_2=44;cs_2=25;re_2=38;ce_2=87;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_3=38;cs_3=87;re_3=34;ce_3=118;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
rs_4=34;cs_4=118;re_4=24;ce_4=178;[path_4, mat_4, ~]=shortpath(re_4,ce_4,rs_4,cs_4,mat);%% calcul de la trajectoire de base  => distance la plus courte entre le point initial et le point final en evitant les point sans dates de colonisation et deja utiliser                  
path=[path_1;path_2;path_3;path_4];path(1:length(path(:,1)))=(1:length(path(:,1)))';mat.flag=mat_1.flag+mat_2.flag+mat_3.flag+mat_4.flag;[row,col]=find(mat.flag>1);mat.flag(row(:),col(:))=1;
coast = load('coast');seaclr=brewermap(15,'RdBu');seaclr1=seaclr(9,:);
figure('Color','w'),ax=worldmap({'Europe','Asia'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_3,cs_3),mat.lon(rs_3,cs_3),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                  
                    
%% calcul de l'entropie de base  => somme des proba entre chaque cellule de grille
np=length(path(:,1));p=proba3(path);inipath=path;Entropie=sum(p);clear p;
stop=0;TMH=20;MHrate=2.7;MHent=[];TMHent=[];stoptime=20000;cpt=1;%0.12 est pas mal
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


figure('Color','w'),ax=worldmap({'Europe','Asia'});gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs_2,cs_2),mat.lon(rs_2,cs_2),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re_1,ce_1),mat.lon(re_1,ce_1),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);

figure('Color','w'),ax=worldmap('World');setm(ax, 'Origin', [0 160 0]);gridm('on');setm(gca,'ParallelLabel','off','MeridianLabel','off');
                    geoshow(coast.lat,coast.long,'DisplayType','polygon','Facecolor','w'),
                    geoshow(flipud(coast.lat),flipud(coast.long),'DisplayType','polygon','FaceColor',seaclr1),
                    hold on,surfm(yg+0.5,xg2-0.5,mat.time);hold on,plotm(mat.lat(rs,cs),mat.lon(rs,cs),'dm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(mat.lat(re,ce),mat.lon(re,ce),'pm','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    hold on,plotm(path(:,4),path(:,5),'-^m','Markersize',6,'Markerfacecolor','m','Markeredgecolor',[0 0 0]);
                    
save -ascii Trajectory_Africa-Beringia_(path#6).txt path;save -ascii Trajectory_Africa-Beringia_(enthropy#6).txt MHent;
save -ascii Trajectory_Africa-Beringia_(Temp#6).txt TMHent;save -ascii Trajectory_Africa-Beringia_(Finalproba#6).txt pfinal;

%%=========================================================================
%% SHORT PATHS
%%=========================================================================
% compilation de toutes les trajectoires initiales avant minimisation
% short path #1 => the shortest to use
clear all, close all, clc,
mat.time=load('Humancolonisation_timing_v4.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs=60;cs=62;re=24;ce=178;[path, mat, ~]=shortpath(re,ce,rs,cs,mat);p=proba3(path);Entropie=sum(p);
save -ascii Short_Trajectory_Africa-Beringia_(path#1).txt path;save -ascii Short_Trajectory_Africa-Beringia_(enthropy#1).txt Entropie;save -ascii Short_Trajectory_Africa-Beringia_(Finalproba#1).txt p;

% short path #2  
clear all, close all, clc,
mat.time=load('Humancolonisation_timing_v4.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=60;cs_1=62;re_1=48;ce_1=87;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);rs_2=48;cs_2=87;re_2=24;ce_2=178;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);
path=[path_1;path_2];p=proba3(path);Entropie=sum(p);
save -ascii Short_Trajectory_Africa-Beringia_(path#2).txt path;save -ascii Short_Trajectory_Africa-Beringia_(enthropy#2).txt Entropie;save -ascii Short_Trajectory_Africa-Beringia_(Finalproba#2).txt p;

% short path #3  
clear all, close all, clc,
mat.time=load('Humancolonisation_timing_v4.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=60;cs_1=62;re_1=44;ce_1=127;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);rs_2=44;cs_2=127;re_2=24;ce_2=178;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);
path=[path_1;path_2];p=proba3(path);Entropie=sum(p);
save -ascii Short_Trajectory_Africa-Beringia_(path#3).txt path;save -ascii Short_Trajectory_Africa-Beringia_(enthropy#3).txt Entropie;save -ascii Short_Trajectory_Africa-Beringia_(Finalproba#3).txt p;

% short path #4  
clear all, close all, clc,
mat.time=load('Humancolonisation_timing_v4.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=60;cs_1=62;re_1=76;ce_1=125;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);rs_2=76;cs_2=125;re_2=24;ce_2=178;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);
path=[path_1;path_2];p=proba3(path);Entropie=sum(p);
save -ascii Short_Trajectory_Africa-Beringia_(path#4).txt path;save -ascii Short_Trajectory_Africa-Beringia_(enthropy#4).txt Entropie;save -ascii Short_Trajectory_Africa-Beringia_(Finalproba#4).txt p;

% short path #5  
clear all, close all, clc,
mat.time=load('Humancolonisation_timing_v4.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=60;cs_1=62;re_1=58;ce_1=114;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);rs_2=58;cs_2=114;re_2=46;ce_2=135;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);
rs_3=46;cs_3=135;re_3=24;ce_3=178;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);path=[path_1;path_2;path_3];p=proba3(path);Entropie=sum(p);
save -ascii Short_Trajectory_Africa-Beringia_(path#5).txt path;save -ascii Short_Trajectory_Africa-Beringia_(enthropy#5).txt Entropie;save -ascii Short_Trajectory_Africa-Beringia_(Finalproba#5).txt p;

% short path #6  
clear all, close all, clc,
mat.time=load('Humancolonisation_timing_v4.txt');mat.time(35:36,25:35)=NaN;[l,c]=size(mat.time);mat.flag=zeros(l,c);
tr=160;xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');mat.lat=repmat(yg,1,c);mat.lon=repmat(xg2,l,1);
rs_1=60;cs_1=62;re_1=44;ce_1=25;[path_1, mat_1, ~]=shortpath(re_1,ce_1,rs_1,cs_1,mat);rs_2=44;cs_2=25;re_2=38;ce_2=87;[path_2, mat_2, ~]=shortpath(re_2,ce_2,rs_2,cs_2,mat);
rs_3=38;cs_3=87;re_3=34;ce_3=118;[path_3, mat_3, ~]=shortpath(re_3,ce_3,rs_3,cs_3,mat);rs_4=34;cs_4=118;re_4=24;ce_4=178;[path_4, mat_4, ~]=shortpath(re_4,ce_4,rs_4,cs_4,mat);path=[path_1;path_2;path_3;path_4];p=proba3(path);Entropie=sum(p);
save -ascii Short_Trajectory_Africa-Beringia_(path#6).txt path;save -ascii Short_Trajectory_Africa-Beringia_(enthropy#6).txt Entropie;save -ascii Short_Trajectory_Africa-Beringia_(Finalproba#6).txt p;





