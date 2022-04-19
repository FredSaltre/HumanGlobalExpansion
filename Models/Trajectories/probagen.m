function p=probagen(path,pw)

np=length(path(:,1));c=5;p=[];%c est une constante de penalite
for i=2:np;
    delta=pw(pw(:,1)==path(i-1,6),pw(1,:)==path(i,6));
    if delta<0;p=[p;(-delta)+c]; else p=[p;c];end;
    clear delta;
end;

end