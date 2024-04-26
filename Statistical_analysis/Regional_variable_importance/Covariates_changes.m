program Covariates_changes(stats)

% This script helps understand if the environmental conditions along certain human migration paths are significantly different from others and whether these differences are statistically robust or could merely be coincidental. 
% It is designed for analyzing environmental data related to human migration paths. It compares various migration trajectories to determine how different environmental factors, such as temperature or distance to rivers, vary between preferred ("Yes" trajectories) and less preferred ("No" trajectories) migration routes
% the same procedure is used for each region as follow:

%1. Data Loading and Preparation
%The script starts by loading a dataset from a text file that contains environmental variables along different migration paths. These paths are distinguished by identifiers that note whether they are preferred paths or not. The script initially filters out any rows from the dataset that contain missing data (NaNs) for the variables of interest.

%2. Variable Selection
%It focuses on analyzing specific environmental variables such as temperature, precipitation, and distance to rivers. The script sets up to use one of these variables (specified by varclim) for the detailed analysis.

%3. Region-Specific Analysis
%The script then narrows down the data to focus specifically on the region "Africa-Beringia." It separates the migration paths into preferred and less preferred based on their identifiers and calculates the median values of the selected environmental variable for these paths.

%4. Comparison of Paths
%It calculates the differences between the median values of the environmental variable for the preferred and less preferred paths. This involves computing pairwise differences, essentially looking at how much more or less favorable the environmental conditions are on preferred paths compared to others.

%5. Statistical Testing
%The script performs a permutation test, which is a statistical method to determine if the observed differences between the preferred and non-preferred paths are statistically significant or if they could be due to random chance. This involves randomly shuffling the path identifiers and recalculating the differences multiple times to create a distribution of differences that could happen by chance.

%6. Visualization and Results
%Finally, it visualizes these differences using histograms and density plots to provide a visual representation of how significant the differences are. It also calculates and records the probability (p-value) that the observed differences could occur by random chance.


%WE DO PAIRWISE COMPARISON OF THE MEDIAN OF TRAJECTORIES
%file format for the statsitical analysis of environmental variable and trajectories.
%input file = New_HumanTrajEnv_(Veglow) = Vegetation rounded down
% COLUMNS ::: [1] = Main Trajectoire => YES = 1 and 1.5 or NO = 2
%         ::: [2] = Region (1 = AfricaBeringia, 2 = AfricaFinland, 3 = AfricaJapan, 4 = AfricaPortugal, 5 = BeringiaMexico, 6 = MexicoBrazil, 7 = MexicoChile)
%         ::: [3] = Trajectoire ID
%         ::: [4-5] = Latitude and Longitude
%         ::: [6] = timing of human arrival
%         ::: [7] = Temperature
%         ::: [8] = Precipitation
%         ::: [9] = Distance to coast
%         ::: [10] = Type of vegetation (1 = forest, 2 = non forest)
%         ::: [11] = Distance to rivers
%         ::: [12] = Ruggedness
% OUTPUT file for each variable and trajectory we have median value except for the type of vegetation = percentage of forest along each trajectory
% by Frederik Saltre - 2/12/2020

%We are just going to look at Temperature, Precipitation, Ditance to rivers and Vegetation


%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% script for: temperature [mat(:,7)], precipitation [mat(:,8)] and distance to rivers [mat(:,11)] = will be tagged by VARCLIM
%% anomalie relative to the prefered trajectories => positive value = Yes trajectories are warmer, wetter and further away from river than the NO trajectories
%% STATS is the output file = each row = 1 region while column are as follow: [,1] = ID region, [,2:4] = 25, 50 and 75 percentile of pairwise, [,5] = pvalue randomisation
clear all, close all, clc, % Clears all variables, closes all figures, and clears the command window

mat=load('New_HumanTrajEnv_(Veglow).txt');varclim = 8; % Load data from a text file and initialize variable for analysis & Index of the variable of interest (e.g., Precipitation)
mat=[mat(:,1:5),mat(:,varclim)];TF=isnan(mat);[row,~]=find(TF==1);mat(row,:)=[];clear TF row;%% Filter out rows with NaN values in the column of interest

% Initialize color schemes for plotting
clr1=brewermap(8,'Pastel2');clr2=brewermap(8,'Set2');clr3=brewermap(8,'Dark2');
% Initialize storage for statistics results
stats=nan(7,5);idstat=0; % 7 regions, 5 stats metrics & Index for stats array

Fig2_DS=[];% Initialize an empty table for storing pairwise distance data

%%=========================================================================
%% STEP #1: AFRICA-BERINGIA
%%=========================================================================
% Processing for AFRICA-BERINGIA region
m1=mat(mat(:,2)==1,:);idstat=idstat+1; %Select data for Africa-Beringia region
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2; % Identify "Yes" trajectories (1 and 1.5 in the dataset)
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[]; % Extract data for "Yes" trajectories and others
idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[]; % Unique trajectory IDs for "Yes" and "No"

% Compute pairwise differences for the specified variable
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
        out=[out;ClimY-ClimN];clear tempN;
    end;clear tempY;  
end;
% Density estimation of pairwise differences
[f_dif1,xi_dif1] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 

% Create a table to store trajectory information
nout=length(out(:,1));
chars = repmat('AFRICA-BERINGIA',nout,1);
Fig2_DS = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Fig2_DS.Trajectory = string(Fig2_DS.Trajectory);% Convert character arrays to strings

clear out idy ym1 nm1 idyes idno;
figure(1),area(xi_dif1,f_dif1,'Facecolor',clr1(1,:),'Edgecolor',clr2(1,:)), % Plot the density distribution

%% Permutation test to calculate p-values
nrandmax=10000; % Number of permutations
resdist=nan(nrandmax,1);% Storage for test statistics
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    % Recalculate pairwise distances after shuffling (similar block of code as above for recalculating 'out' array)
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

% Plot histogram of permutation results and observed statistic
figure(2),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('AFRICA & BERINGIA - PRECIPITATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;

%%=========================================================================
%% STEP #2: AFRICA-FINLAND
%%=========================================================================
%distribution of pairwise difference
m1=mat(mat(:,2)==2,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
        out=[out;ClimY-ClimN];clear tempN;
    end;clear tempY;  
end;[f_dif2,xi_dif2] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 

nout=length(out(:,1));
chars = repmat('AFRICA-FINLAND',nout,1);
Figtemp = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Figtemp.Trajectory = string(Figtemp.Trajectory);
Fig2_DS =[Fig2_DS;Figtemp];clear Figtemp;

clear out idy ym1 nm1 idyes idno;
figure(3),area(xi_dif2,f_dif2,'Facecolor',clr1(2,:),'Edgecolor',clr2(2,:)),



%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(4),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('AFRICA-FINLAND - PRECIPITATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;




%%=========================================================================
%% STEP #3: AFRICA-JAPAN
%%=========================================================================
%distribution of pairwise difference 
m1=mat(mat(:,2)==3,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
        out=[out;ClimY-ClimN];clear tempN;
    end;clear tempY;  
end;[f_dif3,xi_dif3] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 

nout=length(out(:,1));
chars = repmat('AFRICA-JAPAN',nout,1);
Figtemp = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Figtemp.Trajectory = string(Figtemp.Trajectory);
Fig2_DS =[Fig2_DS;Figtemp];clear Figtemp;


clear out idy ym1 nm1 idyes idno;
figure(5),area(xi_dif3,f_dif3,'Facecolor',clr1(3,:),'Edgecolor',clr2(3,:)),



%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(6),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('AFRICA-JAPAN - PRECIPITATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;



%%=========================================================================
%% STEP #4: AFRICA-PORTUGAL
%%=========================================================================
%distribution of pairwise difference 
m1=mat(mat(:,2)==4,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
        out=[out;ClimY-ClimN];clear tempN;
    end;clear tempY;  
end;[f_dif4,xi_dif4] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 

nout=length(out(:,1));
chars = repmat('AFRICA-PORTUGAL',nout,1);
Figtemp = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Figtemp.Trajectory = string(Figtemp.Trajectory);
Fig2_DS =[Fig2_DS;Figtemp];clear Figtemp;

clear out idy ym1 nm1 idyes idno;
figure(7),area(xi_dif4,f_dif4,'Facecolor',clr1(4,:),'Edgecolor',clr2(4,:)),


%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(8),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('AFRICA-PORTUGAL - PRECIPITATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;


%%=========================================================================
%% STEP #5: BERINGIA-MEXICO
%%=========================================================================
%distribution of pairwise difference 
m1=mat(mat(:,2)==5,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
        out=[out;ClimY-ClimN];clear tempN;
    end;clear tempY;  
end;[f_dif5,xi_dif5] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 

nout=length(out(:,1));
chars = repmat('BERINGIA-MEXICO',nout,1);
Figtemp = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Figtemp.Trajectory = string(Figtemp.Trajectory);
Fig2_DS =[Fig2_DS;Figtemp];clear Figtemp;

clear out idy ym1 nm1 idyes idno;
figure(9),area(xi_dif5,f_dif5,'Facecolor',clr1(5,:),'Edgecolor',clr2(5,:)),


%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(10),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('BERINGIA-MEXICO - PRECIPITATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;



%%=========================================================================
%% STEP #6: MEXICO-BRAZIL
%%=========================================================================
%distribution of pairwise difference 
m1=mat(mat(:,2)==6,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
        out=[out;ClimY-ClimN]; %anomalie relative to the prefered trajectories => positive value = Yes trajectories are warmer ; negative value = Yes trajectories are cooler
        clear tempN;
    end;clear tempY;  
end;[f_dif6,xi_dif6] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 

nout=length(out(:,1));
chars = repmat('MEXICO-BRAZIL',nout,1);
Figtemp = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Figtemp.Trajectory = string(Figtemp.Trajectory);
Fig2_DS =[Fig2_DS;Figtemp];clear Figtemp;


clear out idy ym1 nm1 idyes idno;
figure(11),area(xi_dif6,f_dif6,'Facecolor',clr1(6,:),'Edgecolor',clr2(6,:)),


%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(12),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('MEXICO-BRAZIL - PRECIPITATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;




%%=========================================================================
%% STEP #8: MEXICO-CHILE
%%=========================================================================
%distribution of pairwise difference 
m1=mat(mat(:,2)==8,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
        out=[out;ClimY-ClimN]; %anomalie relative to the non prefered trajectories => positive value = Yes trajectories are warmer ; negative value = Yes trajectories are cooler
        clear tempN;
    end;clear tempY;  
end;[f_dif7,xi_dif7] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 

nout=length(out(:,1));
chars = repmat('MEXICO-CHILE',nout,1);
Figtemp = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Figtemp.Trajectory = string(Figtemp.Trajectory);
Fig2_DS =[Fig2_DS;Figtemp];clear Figtemp;

clear out idy ym1 nm1 idyes idno;
figure(13),area(xi_dif7,f_dif7,'Facecolor',clr1(7,:),'Edgecolor',clr2(7,:)),

%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);ClimY=median(tempY(:,end)); 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);ClimN=median(tempN(:,end)); 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(14),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('MEXICO-CHILE - PRECIPITATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;

% Save the table to a CSV file
writetable(Fig2_DS, 'SourceData_Fig2(precipitation).csv');



%%=========================================================================
%% STEP #9: DISPLAY AND SAVING OUTPUTS
%%=========================================================================
%% savings stats results
save -ascii Dist2Rivers_RegionalChanges_(Stats) stats;

%% displaying all figure
size=0.2;basevalue=-10;%basevalue = -10 for all variables but distance to rivers (becomes = -2000)
% scale => Temperature [-5 15]; Precipitation [-0.5 1.5]; Distance to river [-1000 1500]
figure(15),area(xi_dif1,f_dif1,'Facecolor',clr1(1,:),'Edgecolor',clr2(1,:)),
          hold on,area(xi_dif2,f_dif2,'Facecolor',clr1(2,:),'Edgecolor',clr2(2,:)),
          hold on,area(xi_dif3,f_dif3,'Facecolor',clr1(3,:),'Edgecolor',clr2(3,:)),
          hold on,area(xi_dif4,f_dif4,'Facecolor',clr1(4,:),'Edgecolor',clr2(4,:)),
          hold on,area(xi_dif5,f_dif5,'Facecolor',clr1(5,:),'Edgecolor',clr2(5,:)),
          hold on,area(xi_dif6,f_dif6,'Facecolor',clr1(6,:),'Edgecolor',clr2(6,:)),
          hold on,area(xi_dif7,f_dif7,'Facecolor',clr1(7,:),'Edgecolor',clr2(7,:)),
          alpha 0.5,xlim([-0.5 1.5]),axis square,
          hold on,plot(xi_dif1,f_dif1,'Color',clr2(1,:)),          
          hold on,plot(xi_dif2,f_dif2,'Color',clr2(2,:)),
          hold on,plot(xi_dif3,f_dif3,'Color',clr2(3,:)),
          hold on,plot(xi_dif4,f_dif4,'Color',clr2(4,:)),
          hold on,plot(xi_dif5,f_dif5,'Color',clr2(5,:)),
          hold on,plot(xi_dif6,f_dif6,'Color',clr2(6,:)),
          hold on,plot(xi_dif7,f_dif7,'Color',clr2(7,:)),


figure(16),barh(stats(1,1),stats(1,4),size,'basevalue',basevalue,'Facecolor',clr1(1,:),'Edgecolor','none'),hold on,barh(stats(1,1),stats(1,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(1,3),stats(1,1),'|k','Markersize',7,'color',clr3(1,:))
          hold on,barh(stats(2,1),stats(2,4),size,'basevalue',basevalue,'Facecolor',clr1(2,:),'Edgecolor','none'),hold on,barh(stats(2,1),stats(2,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(2,3),stats(2,1),'|k','Markersize',7,'color',clr3(2,:))
          hold on,barh(stats(3,1),stats(3,4),size,'basevalue',basevalue,'Facecolor',clr1(3,:),'Edgecolor','none'),hold on,barh(stats(3,1),stats(3,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(3,3),stats(3,1),'|k','Markersize',7,'color',clr3(3,:))
          hold on,barh(stats(4,1),stats(4,4),size,'basevalue',basevalue,'Facecolor',clr1(4,:),'Edgecolor','none'),hold on,barh(stats(4,1),stats(4,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(4,3),stats(4,1),'|k','Markersize',7,'color',clr3(4,:))
          hold on,barh(stats(5,1),stats(5,4),size,'basevalue',basevalue,'Facecolor',clr1(5,:),'Edgecolor','none'),hold on,barh(stats(5,1),stats(5,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(5,3),stats(5,1),'|k','Markersize',7,'color',clr3(5,:))
          hold on,barh(stats(6,1),stats(6,4),size,'basevalue',basevalue,'Facecolor',clr1(6,:),'Edgecolor','none'),hold on,barh(stats(6,1),stats(6,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(6,3),stats(6,1),'|k','Markersize',7,'color',clr3(6,:))
          hold on,barh(stats(7,1),stats(7,4),size,'basevalue',basevalue,'Facecolor',clr1(7,:),'Edgecolor','none'),hold on,barh(stats(7,1),stats(7,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(7,3),stats(7,1),'|k','Markersize',7,'color',clr3(7,:))
          xlim([-0.5 1.5]),axis square,
          






%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% script for: vegetation [mat(:,10)], calculation on percentage of forest (not the median)= will be tagged by VARCLIM
%% anomalie relative to the prefered trajectories => positive value = Yes trajectories are warmer, wetter and further away from river than the NO trajectories
%% STATS is the output file = each row = 1 region while column are as follow: [,1] = ID region, [,2:4] = 25, 50 and 75 percentile of pairwise, [,5] = pvalue randomisation
%%=========================================================================
clear all, close all, clc,
mat=load('New_HumanTrajEnv_(Veglow).txt');varclim = 10;
mat=[mat(:,1:5),mat(:,varclim)];TF=isnan(mat);[row,~]=find(TF==1);mat(row,:)=[];clear TF row;%idreg=unique(mat(:,2));
clr1=brewermap(8,'Pastel2');clr2=brewermap(8,'Set2');clr3=brewermap(8,'Dark2');stats=nan(7,5);idstat=0;

Fig2_DS=[];

%%=========================================================================
%% STEP #1: AFRICA-BERINGIA
%%=========================================================================
m1=mat(mat(:,2)==1,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
        out=[out;ClimY-ClimN];clear tempN;
    end;clear tempY;  
end;[f_dif1,xi_dif1] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 

nout=length(out(:,1));
chars = repmat('AFRICA-BERINGIA',nout,1);
Fig2_DS = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
% Convert character arrays to string arrays
Fig2_DS.Trajectory = string(Fig2_DS.Trajectory);

clear out idy ym1 nm1 idyes idno;
figure(1),area(xi_dif1,f_dif1,'Facecolor',clr1(1,:),'Edgecolor',clr2(1,:)),


%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(2),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('AFRICA & BERINGIA - VEGETATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;



%%=========================================================================
%% STEP #2: AFRICA-FINLAND
%%=========================================================================
m1=mat(mat(:,2)==2,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
        out=[out;ClimY-ClimN];clear tempN;
    end;clear tempY;  
end;[f_dif2,xi_dif2] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 


nout=length(out(:,1));
chars = repmat('AFRICA-FINLAND',nout,1);
Figtemp = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Figtemp.Trajectory = string(Figtemp.Trajectory);
Fig2_DS =[Fig2_DS;Figtemp];clear Figtemp;


clear out idy ym1 nm1 idyes idno;
figure(3),area(xi_dif2,f_dif2,'Facecolor',clr1(2,:),'Edgecolor',clr2(2,:)),


%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(4),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('AFRICA & FINLAND - VEGETATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;



%%=========================================================================
%% STEP #3: AFRICA-JAPAN
%%=========================================================================
m1=mat(mat(:,2)==3,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
        out=[out;ClimY-ClimN];clear tempN;
    end;clear tempY;   
end;[f_dif3,xi_dif3] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 


nout=length(out(:,1));
chars = repmat('AFRICA-JAPAN',nout,1);
Figtemp = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Figtemp.Trajectory = string(Figtemp.Trajectory);
Fig2_DS =[Fig2_DS;Figtemp];clear Figtemp;

clear out idy ym1 nm1 idyes idno;
figure(5),area(xi_dif3,f_dif3,'Facecolor',clr1(3,:),'Edgecolor',clr2(3,:)),



%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(6),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('AFRICA & JAPAN - VEGETATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;



%%=========================================================================
%% STEP #4: AFRICA-PORTUGAL
%%=========================================================================
m1=mat(mat(:,2)==4,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
        out=[out;ClimY-ClimN]; %anomalie relative to the non prefered trajectories => positive value = Yes trajectories are warmer ; negative value = Yes trajectories are cooler
        clear tempN;
    end;clear tempY; 
end;[f_dif4,xi_dif4] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 

nout=length(out(:,1));
chars = repmat('AFRICA-PORTUGAL',nout,1);
Figtemp = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Figtemp.Trajectory = string(Figtemp.Trajectory);
Fig2_DS =[Fig2_DS;Figtemp];clear Figtemp;

clear out idy ym1 nm1 idyes idno;
figure(7),area(xi_dif4,f_dif4,'Facecolor',clr1(4,:),'Edgecolor',clr2(4,:)),


%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(8),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('AFRICA & PORTUGAL - VEGETATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;


%%=========================================================================
%% STEP #5: BERINGIA-MEXICO
%%=========================================================================
m1=mat(mat(:,2)==5,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
        out=[out;ClimY-ClimN]; %anomalie relative to the non prefered trajectories => positive value = Yes trajectories are warmer ; negative value = Yes trajectories are cooler
        clear tempN;
    end;clear tempY; 
end;[f_dif5,xi_dif5] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 


nout=length(out(:,1));
chars = repmat('BERINGIA-MEXICO',nout,1);
Figtemp = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Figtemp.Trajectory = string(Figtemp.Trajectory);
Fig2_DS =[Fig2_DS;Figtemp];clear Figtemp;

clear out idy ym1 nm1 idyes idno;
figure(9),area(xi_dif5,f_dif5,'Facecolor',clr1(5,:),'Edgecolor',clr2(5,:)),


%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(10),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('BERINGIA & MEXICO - VEGETATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;





%%=========================================================================
%% STEP #6: MEXICO-BRAZIL
%%=========================================================================
m1=mat(mat(:,2)==6,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
        out=[out;ClimY-ClimN]; %anomalie relative to the non prefered trajectories => positive value = Yes trajectories are warmer ; negative value = Yes trajectories are cooler
        clear tempN;
    end;clear tempY;  
end;[f_dif6,xi_dif6] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 

nout=length(out(:,1));
chars = repmat('MEXICO-BRAZIL',nout,1);
Figtemp = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Figtemp.Trajectory = string(Figtemp.Trajectory);
Fig2_DS =[Fig2_DS;Figtemp];clear Figtemp;

clear out idy ym1 nm1 idyes idno;
figure(11),area(xi_dif6,f_dif6,'Facecolor',clr1(6,:),'Edgecolor',clr2(6,:)),



%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(12),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('MEXICO & BRAZIL - VEGETATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;

%%=========================================================================
%% STEP #7: MEXICO-CHILE
%%=========================================================================
m1=mat(mat(:,2)==8,:);idstat=idstat+1;
idy1=find(m1(:,1)==1);idy2=find(m1(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
ym1=m1(idy,:);nm1=m1;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
for i=1:length(idyes);
    tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
    for j=1:length(idno);
        tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
        out=[out;ClimY-ClimN];clear tempN;
    end;clear tempY; 
end;[f_dif7,xi_dif7] = ksdensity(out);stats(idstat,1:4)=[idstat,quantile(out,[0.25 0.5 0.75])]; 

nout=length(out(:,1));
chars = repmat('MEXICO-CHILE',nout,1);
Figtemp = table(chars, out, 'VariableNames', {'Trajectory', 'Pairwise_Distance'});
Figtemp.Trajectory = string(Figtemp.Trajectory);
Fig2_DS =[Fig2_DS;Figtemp];clear Figtemp;


clear out idy ym1 nm1 idyes idno;
figure(13),area(xi_dif7,f_dif7,'Facecolor',clr1(7,:),'Edgecolor',clr2(7,:)),


%randomisation test for pvaluew calculation
nrandmax=10000; % nb of permutation test
resdist=nan(nrandmax,1);%vector to store the output of each permutation stats
f = waitbar(0,'Please wait... permutation in progress');

%% start the permutation loop
for k=1:nrandmax;waitbar(k/nrandmax);
    m2=m1;%copy of original dataset
    idx = randperm(length(m2(:,1))); %% Generate a random permutation of indices rows of the dataset
    m2(:,1) = m2(idx,1);% Use the indices to reshuffle the vector

    %we recalculate the pairwise distance
    idy1=find(m2(:,1)==1);idy2=find(m2(:,1)==1.5);idy=[idy1;idy2];clear idy1 idy2;
    ym1=m2(idy,:);nm1=m2;nm1(idy,:)=[];idyes=unique(ym1(:,3));idno=unique(nm1(:,3));out=[];
    for i=1:length(idyes);
        tempY=ym1(ym1(:,3)==idyes(i),:);countY=find(tempY(:,end)==1);ClimY=(length(countY).*100)/length(tempY(:,1));clear countY; 
        for j=1:length(idno);
            tempN=nm1(nm1(:,3)==idno(j),:);countN=find(tempN(:,end)==1);ClimN=(length(countN).*100)/length(tempN(:,1));clear countN; 
            out=[out;ClimY-ClimN];clear tempN;
        end;clear tempY;  
    end;resdist(k,1)=median(out,'omitnan');
    clear out idy idx ym1 m2 nm1 idyes idno;
end;close(f);
%%end of permutation loop

figure(14),histogram(resdist);hold on, xline(stats(idstat,2),'--r', 'observation', 'LineWidth', 2);
          title('MEXICO & CHILE - VEGETATION', 'FontSize',14);axis square;xlab('median pairwise dist random',13);

%pvalue calculation
idpv=find(abs(resdist)>=abs(stats(idstat,2)));stats(idstat,5)=length(idpv)/nrandmax;
clear m1 idpv resdist;


% Save the table to a CSV file
writetable(Fig2_DS, 'SourceData_Fig2(vegetation).csv');


%%=========================================================================
%% STEP #8: DISPLAY AND SAVING OUTPUTS
%%=========================================================================
%% savings stats results
save -ascii Vegetation_RegionalChanges_(Stats) stats;

%% displaying all figure
size=0.2;basevalue=-200;%basevalue = 200
% scale => Temperature [-10 15]; Precipitation [-1.5 0.5]; Distance to river [-1500 1000]
figure(8),area(xi_dif1,f_dif1,'Facecolor',clr1(1,:),'Edgecolor',clr2(1,:)),
          hold on,area(xi_dif2,f_dif2,'Facecolor',clr1(2,:),'Edgecolor',clr2(2,:)),
          hold on,area(xi_dif3,f_dif3,'Facecolor',clr1(3,:),'Edgecolor',clr2(3,:)),
          hold on,area(xi_dif4,f_dif4,'Facecolor',clr1(4,:),'Edgecolor',clr2(4,:)),
          hold on,area(xi_dif5,f_dif5,'Facecolor',clr1(5,:),'Edgecolor',clr2(5,:)),
          hold on,area(xi_dif6,f_dif6,'Facecolor',clr1(6,:),'Edgecolor',clr2(6,:)),
          hold on,area(xi_dif7,f_dif7,'Facecolor',clr1(7,:),'Edgecolor',clr2(7,:)),
          alpha 0.5,xlim([-100 100]),axis square,
          hold on,plot(xi_dif1,f_dif1,'Color',clr2(1,:)),          
          hold on,plot(xi_dif2,f_dif2,'Color',clr2(2,:)),
          hold on,plot(xi_dif3,f_dif3,'Color',clr2(3,:)),
          hold on,plot(xi_dif4,f_dif4,'Color',clr2(4,:)),
          hold on,plot(xi_dif5,f_dif5,'Color',clr2(5,:)),
          hold on,plot(xi_dif6,f_dif6,'Color',clr2(6,:)),
          hold on,plot(xi_dif7,f_dif7,'Color',clr2(7,:)),


figure(9),barh(stats(1,1),stats(1,4),size,'basevalue',basevalue,'Facecolor',clr1(1,:),'Edgecolor','none'),hold on,barh(stats(1,1),stats(1,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(1,3),stats(1,1),'|k','Markersize',7,'color',clr3(1,:))
          hold on,barh(stats(2,1),stats(2,4),size,'basevalue',basevalue,'Facecolor',clr1(2,:),'Edgecolor','none'),hold on,barh(stats(2,1),stats(2,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(2,3),stats(2,1),'|k','Markersize',7,'color',clr3(2,:))
          hold on,barh(stats(3,1),stats(3,4),size,'basevalue',basevalue,'Facecolor',clr1(3,:),'Edgecolor','none'),hold on,barh(stats(3,1),stats(3,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(3,3),stats(3,1),'|k','Markersize',7,'color',clr3(3,:))
          hold on,barh(stats(4,1),stats(4,4),size,'basevalue',basevalue,'Facecolor',clr1(4,:),'Edgecolor','none'),hold on,barh(stats(4,1),stats(4,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(4,3),stats(4,1),'|k','Markersize',7,'color',clr3(4,:))
          hold on,barh(stats(5,1),stats(5,4),size,'basevalue',basevalue,'Facecolor',clr1(5,:),'Edgecolor','none'),hold on,barh(stats(5,1),stats(5,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(5,3),stats(5,1),'|k','Markersize',7,'color',clr3(5,:))
          hold on,barh(stats(6,1),stats(6,4),size,'basevalue',basevalue,'Facecolor',clr1(6,:),'Edgecolor','none'),hold on,barh(stats(6,1),stats(6,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(6,3),stats(6,1),'|k','Markersize',7,'color',clr3(6,:))
          hold on,barh(stats(7,1),stats(7,4),size,'basevalue',basevalue,'Facecolor',clr1(7,:),'Edgecolor','none'),hold on,barh(stats(7,1),stats(7,2),size,'basevalue',basevalue,'Facecolor','w','Edgecolor','none'),hold on,plot(stats(7,3),stats(7,1),'|k','Markersize',7,'color',clr3(7,:))
          xlim([-100 100]),axis square,



