function [spath,newmat,noval]=shortpath(re,ce,rs,cs,mat)

    mat2=mat;mat2.flag(re,ce)=0;
    path2=[0,rs,cs,mat2.lat(rs,cs),mat2.lon(rs,cs),mat2.time(rs,cs),1,distance(mat2.lat(re,ce),mat2.lon(re,ce),mat2.lat(rs,cs),mat2.lon(rs,cs))];mat2.flag(path2(end,2),path2(end,3))=1;dist=path2(end,end);
    dist=distance(mat.lat(re,ce),mat.lon(re,ce),mat.lat(rs,cs),mat.lon(rs,cs));%calcul de la distance entre les points de depart et arriver
    while dist>0%on essai tant que l'on a pas atteint la destination
        ngh=neighbourg(path2(end,:),mat2);
        lngh=length(ngh(:,1));ngh=[ngh,nan(lngh,1)];for i=1:lngh;ngh(i,end)=distance(mat.lat(re,ce),mat.lon(re,ce),ngh(i,3),ngh(i,4));end;
        TFngh=isnan(ngh(:,5));if ~isempty(TFngh>0);ngh(TFngh==1,:)=[];end;TFngh2=find(ngh(:,6)==1);if ~isempty(TFngh2>0);ngh(TFngh2,:)=[];end;

        nghflag=0;if isempty(ngh==1);nghflag=1;break;end;
        id=find(ngh(:,end)==min(ngh(:,end)));path2=[path2;path2(end,1)+1,ngh(id(1),:)];
        mat2.flag(path2(end,2),path2(end,3))=1;path2(end,7)=1;%on marque la cellule de grille comme utilisee
        dist=path2(end,end);
        clear ngh id TFngh;   
    end;
    spath=path2;
    newmat=mat2;
    noval=nghflag;       
end