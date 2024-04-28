program Environmental_variables
% developed by Frederik Saltre 26/04/2024
% The script is designed to compute the distance of each grid cell in a specified area to the nearest coast and the nearest river, as well as calculate a ruggedness index for each grid cell
% It includes data preparation, intensive calculations for distance and ruggedness, and sophisticated visualizations with geographical mappings.


%%=========================================================================
%% RETURN A RASTER WITH THE DISTANCE OF EACH GRIDCELL TO THE NEAREST COAST
%%=========================================================================
%Part 1: Distance to the Nearest Coast
%Data Loading: The script starts by loading a matrix from a text file that likely contains geographic or environmental data.
%Preparation: It prepares two sets of data: one for non-missing data points and one specifically for coastal data points.
%Distance Calculation: For each non-missing data point, it calculates the distance to the nearest coast using geographic coordinates. It uses a function to compute the distance between points on the Earth's surface and selects the minimum distance for each point.
%Result Storage: The calculated distances are stored in a matrix that matches the original grid's dimensions and then saved to a text file.

clear all, close all, clc, %% Clear workspace, close all figures, and clear command window
%% Load geographic raster data and prepare for processing
mat=load('maskedHumantimingraster_v2.txt');% Load raster data for human timing
[r,c]=size(mat);tr=160;% Get the dimensions of the matrix + Threshold value (possibly for longitude adjustment)
xg=-179.5:1:179.5;yg=-89.5:89.5;xg2=[xg(:,tr:end),xg(:,1:tr-1)];yg=flipud(yg');% Define longitude and latitude range + Adjust longitudes by threshold + Flip latitude array to match geographic conventions
%% Prepare matrix for non-missing and coastal data
TF=isnan(mat);[row,col]=find(TF==0);nr=length(row);out=[row,col,nan(nr,3)];for i=1:nr;out(i,3:4)=[yg(row(i)),xg2(col(i))];end;% Identify NaNs in the matrix + Find indices of non-NaN entries + Number of non-NaN entries + Initialize output matrix with NaNs for additional data + Populate latitudes and longitudes
%% Prepare coastal data matrix
[row2,col2]=find(TF==1);nr2=length(row2);coastal=[row2,col2,nan(nr2,2)];for i=1:nr2;coastal(i,3:4)=[yg(row2(i)),xg2(col2(i))];end; %similar process as coastal data
%% Calculate distance from each grid cell to the nearest coast
h = waitbar(0,'Please wait...');
for i=1:nr;waitbar(i/nr)
    [vc,~]=distance(out(i,3),out(i,4),coastal(:,3),coastal(:,4));% Calculate distances to coast
    id=find(vc==min(vc));out(i,5)=deg2km(vc(id(1)));% Find the closest coast + Convert closest distance to kilometers
    clear vc id;
end;close(h);
%% Save the distance data to a raster matrix and output to a text file
res=nan(r,c);for i=1:nr;res(out(i,1),out(i,2))=out(i,5);end;% Initialize the result matrix with NaNs + Populate the result matrix with calculated distances
save -ascii NewDistance2coast.txt res;

%% Visualization setup
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
%Part 2: Distance to the Nearest River
%Data and Preparation: Similar to the coast distance calculation, it loads main river data and prepares matrices for calculation.
%Distance Calculation: For each grid cell, it calculates the distance to the nearest river in a similar manner to the coast distance calculation.
%Visualization and Saving: The results are visualized on a world map using MATLAB's mapping toolboxes, and the distance data is saved to a text file.

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
%Part 3: Ruggedness Index Calculation
%Data Loading: This section loads elevation data from a NetCDF file, which contains latitude, longitude, and elevation data.
%Preparation: Adjusts the longitude data to match the grid orientation used in the previous calculations.
%Ruggedness Calculation: For each grid cell, it calculates the ruggedness based on elevation differences with its immediate neighbors. This involves computing the square of the difference in elevation between a cell and each of its eight neighbors, summing these squares, and then taking the square root of the total.
%Final Adjustments and Saving: Excludes areas covered by ice sheets or water bodies from the ruggedness calculation and saves the ruggedness data to a text file.

% This script calculates a ruggedness index based on elevation differences
% between each grid cell and its neighbors.

% Clear workspace, close all figures, and clear command window
clear all, close all, clc,

% Load elevation data from a NetCDF file
mat.alt = ncread('elev.1-deg.nc', 'data');  % Elevation data
mat.lat = ncread('elev.1-deg.nc', 'lat');   % Latitude data
mat.lon = ncread('elev.1-deg.nc', 'lon');   % Longitude data

% Load timing data which may be used to mask certain areas
mat.time = load('maskedHumantimingraster_v2.txt');

% Define the grid adjustment threshold and grid definitions
tr = 160;  % Threshold for adjusting grid alignment
xg = -179.5:1:179.5;  % Longitude grid definition from -179.5 to 179.5
yg = -89.5:89.5;       % Latitude grid definition from -89.5 to 89.5
xg2 = [xg(:, tr:end), xg(:, 1:tr-1)];  % Adjusted longitude grid
yg = flipud(yg');  % Flip latitude grid to match geographic orientation

% Transpose and adjust longitude data to match grid orientation
mat.alt = mat.alt';
mat.lon = [0.5:179.5, -179.5:1:-0.5];  % Adjust longitude data
id = find(mat.lon == xg2(1));
mat.alt2 = [mat.alt(:, id:end), mat.alt(:, 1:id-1)];  % Adjust altitude data grid

% Initialize ruggedness matrix and variables for processing
nr = 180;  % Number of latitude rows
nc = 360;  % Number of longitude columns
cpt = 0;   % Counter for progress tracking
mat.rug = nan(nr, nc);  % Initialize ruggedness matrix with NaNs
h = waitbar(0, 'Please wait...');  % Initialize wait bar for user feedback

% Calculate ruggedness for each grid cell
for i = 2:nr-1
    for j = 2:nc-1
        ngh = nan(9,1);  % Initialize neighbor differences matrix
        cpt = cpt + 1;   % Increment progress counter
        waitbar(cpt / (nr * nc), h);  % Update wait bar
        
        % Calculate squared differences of elevation with eight neighbors
        ngh(1,1) = (mat.alt2(i-1, j-1) - mat.alt2(i, j))^2;
        ngh(2,1) = (mat.alt2(i-1, j) - mat.alt2(i, j))^2;
        ngh(3,1) = (mat.alt2(i-1, j+1) - mat.alt2(i, j))^2;
        ngh(4,1) = (mat.alt2(i, j-1) - mat.alt2(i, j))^2;
        ngh(5,1) = (mat.alt2(i, j) - mat.alt2(i, j))^2;  % Self difference (always 0)
        ngh(6,1) = (mat.alt2(i, j+1) - mat.alt2(i, j))^2;
        ngh(7,1) = (mat.alt2(i+1, j-1) - mat.alt2(i, j))^2;
        ngh(8,1) = (mat.alt2(i+1, j) - mat.alt2(i, j))^2;
        ngh(9,1) = (mat.alt2(i+1, j+1) - mat.alt2(i, j))^2;
        
        % Calculate and store the ruggedness index as the square root of the sum of squared differences
        mat.rug(i, j) = sqrt(sum(ngh));
        clear ngh;
    end;
end;
close(h);  % Close the wait bar

% Mask ruggedness data for areas with NaNs in timing data
TF = isnan(mat.time);  % Identify NaNs in timing data
[row, col] = find(TF == 1);  % Find indices of NaN entries
for i = 1:length(row)
    mat.rug(row(i), col(i)) = NaN;  % Set ruggedness to NaN where timing data is missing
end;
out=mat.rug;save -ascii NewRuggedness.txt out;


