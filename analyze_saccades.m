% This functions analyzes the saccades from viewing Hollywood clip
% conditions

load 'dbclip.mat';
%load 'saccadesWithTrial.mat';

load allSaccades

load baselineSubjectIndex
load subjectIndex
%load subjects_age
close all



% code written after pasting columns from excel with subject id and sex
% % 
% for i=1:length(subjectList)
%     for j=1:length(a)
%         if strcmp(char(subjectList{i,1}), a{j,1})
%             switch (a{j,2})
%                 case 'Male' 
%                     subjectList{i,3}= 1;
%                 case 'Female' 
%                     subjectList{i,3}=0;
%             end
%                         
%             switch char(a{j,3})
%                
%                 case 'Some high school'
%                     subjectList{i,4} = 1;
%                     t(1) = 1;
%                 case 'Some college, no degree'
%                      subjectList{i,4} = 3;
%                      t(3) = 1;
%                 case 'High school diploma'
%                     subjectList{i,4} = 2;
%                     t(2) = 1;
%                 case 'Bachelor''s degree'
%                     subjectList{i,4}= 4;
%                     t(4) = 1;
%                 case 'Associate degree'
%                     subjectList{i,4} = 5;
%                     t(5) = 1;
%                 case 'Master''s degree'
%                     subjectList{i,4}= 6;
%                     t(6) = 1;
%                 case 'professional'
%                     subjectList{i,4} = 7;
%                     t(7) = 1;
%                 case 'Doctoral degree'
%                     subjectList{i,4} = 8;
%                     t(8) = 1;
%                 case 'Primary degree'
%                     subjectList{i,4} = 1;
%                     
%                 otherwise
%                     subjectList{i,4} = nan;
%                     
%             end
%             
%         end
%     end
% end


% Include the age and sex (0 female, 1 male) in the saccade information
j = 1;
for i=10000:length(saccadeIndex)
    
    saccadeIndex(j).age = subjectList{ saccadeIndex(i).subNumber, 2};
    saccadeIndex(j).sex = subjectList{ saccadeIndex(i).subNumber, 3};
    saccadeIndex(j).education = subjectList{ saccadeIndex(i).subNumber, 4};
    magnitude(j) = saccadeIndex(i).length;
    j = j +1;
end

clear medcvn err magnitude

% Extract information about clip, trial, and saccade features so we can
% analyze them and reassign them properly later

load('saccadeTV1Degree_plus16msbefore');
%saccadeIndex = sel_saccades;
magnitude = [];
peakvel = [];
medcvn = [];
err = [];
duration = [];
k=1;
saccadeIndex = saccadeIndexCoilPlus16before;

j = 1;
for k=1:length(saccadeIndex)
    
    %if length(saccadeIndex(j).elapsedtime)>15
        % if ~isempty(saccadeIndex(j).age)
        % age(j) = saccadeIndex(j).age;
        
        % subjects(j) = saccadeIndex(j).subNumber;
        % trials(j) = saccadeIndex(j).trialNumber;
        medcvn(k) = saccadeIndex(k).curvature;
        err(k) = saccadeIndex(k).err;
        direction(k) = saccadeIndex(k).direction;
        magnitude(k) = saccadeIndex(k).length;
        peakvel(k) = saccadeIndex(k).pkvelocity;
        duration(k) = saccadeIndex(k).duration;
        %peakvel(k) = max(saccadeIndex(j).velocity.velocity);
        %j = j+1;
       
   % end
    %     velocity(j) = max(saccadeIndex(j).velocity.velocity);
    
    % end
end


scatter(magnitude, duration, 1);
scatter(magnitude, peakvel, 1);
 S.Axes	= 'Linear';
            S.Minimum_X			= 0;
            S.Maximum_X			= 40;
            S.Minimum_Y			= 0;
            S.Maximum_Y			= 1.5;
            S.Number_of_Colormap_Bins = 600;
            S.Show_textBox      = 1;
            S.Font_Size         = 24;
            S.Marker_Size       = 0.1;
            S.Marker_Line_Width	= 2;
            S.Fit_Line_Width   	= 2;
            S.Restrict_To_X_Range    = 0;
            S.Log_Transform_Data    = 0;
            S.Type_Of_Points  = 'Colormap';
            
plot_mainsequence(magnitude,peakvel,S);

mags= histc(magnitude,0:1:40)./ sum (histc(magnitude,0:1:40));

subplot(1,2,1);
bar(0:1:40,mags,'b','LineWidth', 2); hold on
%hist(fixdurations,[0:100:2500]);
set(gca,'FontSize',18,'FontName','Arial');
xlabel('Saccade magnitude (deg)','FontSize',21,'FontName','Arial');
ylabel('Relative frequency','FontSize',21,'FontName','Arial');
xlim([0 41]);
box off;
subplot(1,2,2);
  axStill = axes();
  h(1) = plot_polardist(axStill, direction);
         


X = [magnitude' peakvel'];
hist3(X, [100 100]);
xlabel('Saccade magnitude (deg)'); ylabel('Saccade peak velocity (deg/ms)');
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
S = findobj('type','surf');
ZD = get(S,'zdata');
ZD(~ZD) = .1;
set(S,'zdata',ZD);
set(gca,'zscale', 'log');

   mYX = [magnitude', peakvel'];
   vYEdge = linspace(0,3,55);
   vXEdge = linspace(0,40,55);
   figure;
   mHist2d = hist2d(mYX,vXEdge,vYEdge); 
   Plot2dHist(mHist2d, vXEdge, vYEdge, 'X', 'Y', '2D Histogram'); colorbar
   xlabel({'Saccade Magnitude (deg)'},'FontSize',18,'FontName','Arial');
   ylabel({'Saccade Peak Velocity (deg/ms)'},'FontSize',18,'FontName','Arial');
   set(gca, 'Xscale',  'log');
   set(gca, 'Yscale',  'log');
 
   
   mYX = [log10(magnitude)', log10(peakvel)'];
   vYEdge = linspace(-1.4,0.3,17);
   vXEdge = linspace(0,1.6,17);
   figure;
   mHist2d = hist2d(mYX,vXEdge,vYEdge);    
  % pcolor (vXEdge(1:end-1), vYEdge(1:end-1), mHist2d);
   
   Plot2dHist(mHist2d, vXEdge, vYEdge, 'X', 'Y', '2D Histogram'); colorbar
   xlabel({'Saccade Magnitude (deg)'},'FontSize',18,'FontName','Arial');
   ylabel({'Saccade Peak Velocity (deg/ms)'},'FontSize',18,'FontName','Arial');
   
   set(gca, 'Xscale',  'log')
   
   
   
   mYX = [log10(magnitude)', log10(peakvel)'];    
   vXEdge = logspace(-2,2,50);
   vYEdge = logspace(-2,1.3,50);
   figure;
   mHist2d = hist2d(mYX,vXEdge,vYEdge); 
   Plot2dHist(mHist2d, vXEdge, vYEdge, 'X', 'Y', '2D Histogram'); colorbar
   
   pcolor( vXEdge, vYEdge,(log10(mHist2d)));
      
   xlabel({'Saccade Magnitude (deg)'},'FontSize',18,'FontName','Arial');
   ylabel({'Saccade Peak Velocity (deg/ms)'},'FontSize',18,'FontName','Arial');
   
% Code written to add some missing velocity information
% for sacind=66787:69272
%    
%     XY=saccadeIndex(sacind).saccade;
%     vel{sacind}=sqrt(sum(diff(XY).^2,2));
%     [velfit{sacind} params{sacind} err(sacind)]=fitVelocity(vel{sacind},4);
%     saccadeIndex(sacind).velocity.velocity=vel{sacind};
%     saccadeIndex(sacind).velocity.velocityFit=velfit{sacind};
%     saccadeIndex(sacind).velocity.params=params{sacind};
%     saccadeIndex(sacind).velocity.err=err(sacind);
% end


% code to calculate the number of total seconds from all subjects

% mysaccades2 is calculated later on , basically converting saccadeIndex
% into a matrix

[unique_subjects order] = unique(mysaccades2(:,3));
trial_count = 0;
trials = [];
for s=1:length(unique_subjects)
    [trials x] = unique(mysaccades2(mysaccades2(:,3) == unique_subjects(s),2));
    trial_count = trial_count + length(trials);    
    
end
 disp(trial_count * 30);

 
 % Study based on subjects
for j=1:length(saccadeIndex)
subj(j) = saccadeIndex(j).subNumber;
clip(j) = saccadeIndex(j).videoNumber;
end

subjects = unique(subj);

% for j=1:length(saccadeIndex)
%     trials(j) = saccadeIndex(j).trialNumber;
% end

% Useful code to assign order 1-40 based on the order subjects watched the
% clip. This will be used to assess fatigue level
[unique_subjects order] = unique(subjects);
my_order= [];
for s=1:length(unique_subjects)
 [unique_trials order] = unique(trials(subjects == unique_subjects(s)));
 %a has the trial number, b has the order (1-40)
 [my_trials my_orders ] = sort(unique_trials)
  for ii=1:length(my_trials)
     my_order(find(trials==my_trials(ii)))=my_orders(ii);      
  end    
end

%saccades after 69272 belongs to defocus study where order is 1-20 and has
%been already included in trialNumber (from code addingDefocusSubjects)
for j=1:length(saccadeIndex)
    
%     if j>69272        
%         saccadeIndex(j).order = saccadeIndex(j).trialNumber;
%     else
       saccadeIndex(j).order = my_order(j);        
%     end
end

%Eliminating some possible noise
magnitude(magnitude>30)=30;
% figure;
% scatter(magnitude(1:66786), velocity(1:66786));
% xlim([0,30]);


% Useful to define the age groups
unique_age = age(order);
y = quantile(unique_age,2);

figure; hist(unique_age,[20:5:90], 'FaceColor',[0.50 0.50 0.50],'EdgeColor',[0.80 0.80 0.80]);
xlim([20,90]);
xlabel('Age','FontSize',24,'FontName','Arial');
% Create ylabel
ylabel('N','FontSize',24,'FontName','Arial');
box off
line([y(1) y(1)], [0 9],'LineStyle','--');
line([y(2) y(2)], [0 9],'LineStyle','--');
% Create title
title({'Age distribution'},'FontSize',24,'FontName','Arial');

% 2D Histogram

% For main sequence
figure
  mYX = [magnitude', peakvel'];
  ndhist(mYX, 'log', 'bins', 100, 'axis', [0 40 0 1.5] , 'normalizeR');  
  scatter(magnitude, peakvel,0.5);
  set(gca, 'xlim', [0 15], 'ylim', [0 1]);
  
  
%    vXEdge = linspace(0,20,100);
%    vYEdge = linspace(0,50,100);
%    mHist2d3 = hist2d(mYX,vYEdge,vXEdge);
%      Plot2dHist(mHist2d3, vXEdge, vYEdge, 'X', 'Y', '2D Histogram'); 
   
  
   mYX = [medcvn', err'];
    ndhist(mYX, 'log', 'bins', 100, 'axis', [0 0.25 0 0.01] , 'normalizeR');
   vXEdge = linspace(0,0.25,100)
   vYEdge = linspace(0,0.01,100);
   mHist2d = hist2d(mYX,vYEdge,vXEdge);
  
   Plot2dHist(mHist2d, vXEdge, vYEdge, 'X', 'Y', '2D Histogram'); colorbar
   xlabel({'Saccade curvature'},'FontSize',18,'FontName','Arial');
   ylabel({'Saccade velocity error'},'FontSize',18,'FontName','Arial');
   
   
   % Plotting the inset 
   figure1 = figure;
   mYX = [magnitude', medcvn'];
   vXEdge = linspace(0,0.25,100);
   vYEdge = linspace(0,30,100);
   mHist2d2 = hist2d(log10(mYX),log10(vXEdge),log10(vYEdge)); 
 
   
   ndhist(mYX, 'log', 'bins', 10, 'axis', [0 30 0 0.25] )%, 'normalizeR');
   %Plot2dHist(mHist2d2, vYEdge, vXEdge, 'X', 'Y', '2D Histogram'); 
  
%    Contours=[1 100 200 400];
%    colorbar('YTick',log10(Contours),'YTickLabel',Contours);
%    colormap(jet);
%    caxis(log10([Contours(1) Contours(length(Contours))]));
%    


xlabel({'Saccade magnitude (deg)'},'FontSize',18,'FontName','Arial');
ylabel({'Saccade peak velocity (deg/s)'},'FontSize',18,'FontName','Arial');
ylabel({'Saccade curvature'},'FontSize',18,'FontName','Arial');
set(gca,'Parent',figure1,'FontSize',16,'FontName','Arial');


figure2 = figure;
mYX = [medcvn', magnitude'];
vXEdge = linspace(0,0.25,100);
vYEdge = linspace(5,30,100);
mHist2d2 = hist2d(mYX,vXEdge,vYEdge);
Plot2dHist(mHist2d2, vYEdge, vXEdge, 'X', 'Y', '2D Histogram');
xlabel({'Saccade magnitude (deg)'},'FontSize',18,'FontName','Arial');
ylabel({'Saccade curvature'},'FontSize',18,'FontName','Arial');
set(gca,'Parent',figure2,'FontSize',16,'FontName','Arial');
   
   
   
   mYX = [magnitude', err'];
   ndhist(mYX, 'log', 'bins', 100, 'axis', [0 30 0 0.01] );%,'normalizeR');
   vXEdge = linspace(0,10,100);
   vYEdge = linspace(0,0.01,100);
   mHist2d3 = hist2d(mYX,vYEdge,vXEdge);
  
   Plot2dHist(mHist2d3, vXEdge, vYEdge, 'X', 'Y', '2D Histogram'); 
   xlabel({'Saccade magnitude'},'FontSize',18,'FontName','Arial');
   ylabel({'Saccade velocity error'},'FontSize',18,'FontName','Arial');
   xlim([0 30]);
   ylim([0 0.01]);
      
   mesh(mHist2d, mHist2d2, mHist2d3)
   xlabel({'Curvature/Error'},'FontSize',16,'FontName','Arial');
   ylabel({'Curvature/Magnitude'},'FontSize',16,'FontName','Arial');
   zlabel({'Magnitude/Error'},'FontSize',16,'FontName','Arial');
   
   % 2D Histogram for landing positions
%    for i=1:length(sac)
%        
%        x=sac{i}(:,2)-sac{i}(1,2);
%        y=sac{i}(:,1)-sac{i}(1,1);       
%        xend(i) = x(end);
%        yend(i) = y(end);
%        
%    end
    for i=1:length(saccadeIndex)
       
       x=saccadeIndex(i).saccade(:,2)-saccadeIndex(i).saccade(1,2);
       y=saccadeIndex(i).saccade(:,1)-saccadeIndex(i).saccade(1,1);       
       xend(i) = x(end);
       yend(i) = y(end);
       
   end
   mYX = [xend', yend'];
   
   mYX(mYX(:,1)>30 | mYX(:,2)>30, :) = [];
   mYX(mYX(:,1)<-30 | mYX(:,2)<-30, :) = [];
   
   ndhist(mYX, 'log', 'bins', 4, 'axis', [-30 30 -30 30] , 'normalizeR');
   
   vXEdge = linspace(-5,5,501);
   vYEdge = linspace(-5,5,501);
   mHist2d3 = hist2d(mYX,vYEdge,vXEdge);
   
   Plot2dHist(mHist2d3, vXEdge, vYEdge, 'X', 'Y', '2D Histogram');
   xlabel({'Horizontal landing point (deg)'},'FontSize',18,'FontName','Arial');
   ylabel({'Vertical landing point (deg)'},'FontSize',18,'FontName','Arial');
     xlim([-30 30]);
   ylim([-30 30]);
   
   
%% This block is pre-processing to modify saccadeIndex into a matlab matrix that
% can be exported into XLS / CSV format to be loaded in Stata
for j=1:length(saccadeIndex)
    
    if isempty(saccadeIndex(j).trialNumber)
        saccadeIndex(j).trialNumber = 0;
    end
    if cuts(saccadeIndex(j).videoNumber) < 5
        saccadeIndex(j).cuts = 0 ;
    else if cuts(saccadeIndex(j).videoNumber) > 5
        saccadeIndex(j).cuts =2 ;
        else 
            saccadeIndex(j).cuts = 1;
        end
    end
    if light05(saccadeIndex(j).videoNumber) < 2
        saccadeIndex(j).light = 0;
    else if light05(saccadeIndex(j).videoNumber) > 3
        saccadeIndex(j).light = 2;
          else 
            saccadeIndex(j).light = 1;
        end
    end
    if audInfo(saccadeIndex(j).videoNumber) < 2
        saccadeIndex(j).audInfo =  0;
    else if audInfo(saccadeIndex(j).videoNumber) > 3
        saccadeIndex(j).audInfo =  2;
         else 
            saccadeIndex(j).audInfo = 1;
        end
    end
    if faces05(saccadeIndex(j).videoNumber) < 2
        saccadeIndex(j).face =  0;
    else if faces05(saccadeIndex(j).videoNumber) > 3
        saccadeIndex(j).face =  2;
          else 
            saccadeIndex(j).face = 1;
        end
    end
    if  strcmp(environment(saccadeIndex(j).videoNumber),'indoor')
        saccadeIndex(j).envi = 0;
    else
        saccadeIndex(j).envi = 1;
    end
    
    % +/- 5 degrees for horizontal . Reference: 
  %  Duke Elder, S. & Wybar, K. (1973). Ocular Motility and Stabismus. (3rd ed.). System of Ophthalmology  (Vol. VI). London: Henry Kimpton.

    if (cos(deg2rad(saccadeIndex(j).direction)) > 0.99) || (cos(deg2rad(saccadeIndex(j).direction)) < -0.99)
        saccadeIndex(j).horver = 1; % horizontal
    else  if  (cos(deg2rad(saccadeIndex(j).direction)) > - 0.3827) && (cos(deg2rad(saccadeIndex(j).direction)) < 0.3827)
            saccadeIndex(j).horver = 3; % vertical
        else
            saccadeIndex(j).horver = 2; % oblique
        end
    end
    
   
    if nature05(saccadeIndex(j).videoNumber) <2
        saccadeIndex(j).nature = 0;
    else if nature05(saccadeIndex(j).videoNumber) >3
        saccadeIndex(j).nature = 2;
        else
             saccadeIndex(j).nature = 1;
        end
    end
   saccadeIndex(j).err = saccadeIndex(j).velocity.err;
   
   if saccadeIndex(j).velocity.velocity(end)>0.3
       saccadeIndex(j).odd = 1;
   else
       saccadeIndex(j).odd = 0;
   end
end

% Report number odd saccades (velocity end > 0.3)
for i=1:length(saccadeIndex)
   my_odd(i) = saccadeIndex(i).odd;    
end
disp(sum(my_odd)/length(my_odd));


saccadeIndex = rmfield(saccadeIndex,'saccade');
saccadeIndex = rmfield(saccadeIndex,'velocity');
saccadeIndex = rmfield(saccadeIndex,'elapsedtime');


for i=1:length(mysaccades2)
    if isnan(mysaccades2(i,9))
        mysaccades2(i,9) = 1;
    end
end

% This code allows our structs with fields to be stored in a excel file
% with just numeric columns
names = fieldnames(saccadePlusCuts);
mysaccades=struct2cell(saccadePlusCuts);
% mysaccades2=cell2mat(mysaccades);

mysaccades = cell2mat(saccades);

mysaccades2 = [];


for i=1:size(mysaccades,1)
    for j=1:size(mysaccades,3)
        saccades(j,i) = mysaccades(i,1,j);
    end
end

for i=1:size(saccades,1)
   
    for j=1:size(saccades,2)
   
        if isempty(cell2mat(saccades(i,j))) || isnan(cell2mat(saccades(i,j)))
            
            mysaccades2(i,j) = nan;
            
        else
           mysaccades2(i,j) = cell2mat(saccades(i,j));
        end
    end    
end


% dlmwrite('saccadeIndex.csv',mysaccades2,'Delimiter',',')
names = names';
xlswrite('saccadePlusCuts.xls',names, 1, 'A1');
xlswrite('SaccadePlusCuts.xls',mysaccades2, 1, 'A1');


%% generic GLOBAL histograms
figure1= figure;
subplot(3,1,1,'Parent',figure1,'FontSize',16,'FontName','Arial');
plot(0:0.01:0.3,histc(medcvn,0:0.01:0.3),'k','LineWidth', 2);
hold on
%plot(0:0.01:0.3,histc(medcvn2,0:0.01:0.3),'r','LineWidth', 2);
xlabel({'Saccade curvature'},'FontSize',18,'FontName','Arial');
ylabel({'N'},'FontSize',18,'FontName','Arial');
xlim([0 0.3]);
box off
subplot(3,1,2,'Parent',figure1,'FontSize',16,'FontName','Arial');
plot(0:0.0005:0.01,histc(err,0:0.0005:0.01),'k','LineWidth', 2);
hold on
%plot(0:0.0005:0.01,histc(err2,0:0.0005:0.01),'r','LineWidth', 2);
ylabel({'N'},'FontSize',18,'FontName','Arial');
xlabel({'Saccade velocity error'},'FontSize',18,'FontName','Arial');
box off
subplot(3,1,3,'Parent',figure1,'FontSize',16,'FontName','Arial');
plot(0:0.5:30,histc(magnitude,0:0.5:30),'k','LineWidth', 2);
hold on
%plot(0:0.5:30,histc(magnitude2,0:0.5:30),'r','LineWidth', 2);
ylabel({'N'},'FontSize',18,'FontName','Arial');
xlabel({'Saccade magnitude'},'FontSize',18,'FontName','Arial');
box off

plot_histograms_odd_saccades(medcvn, medcvn2, err, err2, magnitude, magnitude2, 0, 'Study');


print -dpng histograms

% for j=1:length(saccadeIndex)
% elapsed(j) = saccadeIndex(j).elapsedtime(end);
% end

% Study based on subjects
for j=1:length(saccadeIndex)
subj(j) = saccadeIndex(j).subNumber;
clip(j) = saccadeIndex(j).videoNumber;
end

subjects = unique(subj);
sub_sacc_cvn = [];
sub_sacc_err = [];
sub_sacc_mag = [];
cuts_sacc = [];

for s=1:length(subjects)
   sub_sacc_cvn{s} = medcvn(find(subj==subjects(s)));
   sub_sacc_err{s} = err(find(subj==subjects(s)));
   sub_sacc_mag{s} = magnitude(find(subj==subjects(s)));
end
 
subplot(2,2,1);
for s=1:round(length(subjects)/4)
plot(0:0.01:0.3,histc(sub_sacc_cvn{s},0:0.01:0.3),'LineWidth',2,'Color',[randi([1 100],1,1)/100, randi([1 100],1,1)/100, randi([1 100],1,1)/100]); hold on
end
xlabel({'Saccade curvature'},'FontSize',14,'FontName','Arial');
ylabel({'N'},'FontSize',14,'FontName','Arial');
box off;
set(gca,'XLim',[0 0.25]);
subplot(2,2,2);
for s=round(round(length(subjects)/4)):round(length(subjects)/2)
plot(0:0.01:0.3,histc(sub_sacc_cvn{s},0:0.01:0.3),'LineWidth',2,'Color',[randi([1 100],1,1)/100, randi([1 100],1,1)/100, randi([1 100],1,1)/100]); hold on
end
xlabel({'Saccade curvature'},'FontSize',14,'FontName','Arial');
ylabel({'N'},'FontSize',14,'FontName','Arial');
box off;
set(gca,'XLim',[0 0.25]);
subplot(2,2,3);
for s=round(length(subjects)/2):round(3*length(subjects)/4)
plot(0:0.01:0.3,histc(sub_sacc_cvn{s},0:0.01:0.3),'LineWidth',2,'Color',[randi([1 100],1,1)/100, randi([1 100],1,1)/100, randi([1 100],1,1)/100]); hold on
end
xlabel({'Saccade curvature'},'FontSize',14,'FontName','Arial');
ylabel({'N'},'FontSize',14,'FontName','Arial');
box off;
set(gca,'XLim',[0 0.25]);
subplot(2,2,4);
for s=round(3*length(subjects)/4):length(subjects)
plot(0:0.01:0.3,histc(sub_sacc_cvn{s},0:0.01:0.3),'LineWidth',2,'Color',[randi([1 100],1,1)/100, randi([1 100],1,1)/100, randi([1 100],1,1)/100]); hold on
end
xlabel({'Saccade curvature'},'FontSize',14,'FontName','Arial');
ylabel({'N'},'FontSize',14,'FontName','Arial');
box off;
set(gca,'XLim',[0 0.25]);

print -dpng curvature_subjects


subplot(2,1,1);
for s=1:round(length(subjects)/2)
plot(0:0.0005:0.01,histc(sub_sacc_err{s},0:0.0005:0.01),'LineWidth',2,'Color',[randi([1 100],1,1)/100, randi([1 100],1,1)/100, randi([1 100],1,1)/100] ); hold on;
end
ylabel({'N'},'FontSize',14,'FontName','Arial');
xlabel({'Saccade velocity error'},'FontSize',14,'FontName','Arial');
box off
subplot(2,1,2);
for s=round(length(subjects)/2):length(subjects)
plot(0:0.0005:0.01,histc(sub_sacc_err{s},0:0.0005:0.01),'LineWidth',2,'Color',[randi([1 100],1,1)/100, randi([1 100],1,1)/100, randi([1 100],1,1)/100] ); hold on;
end
ylabel({'N'},'FontSize',14,'FontName','Arial');
xlabel({'Saccade velocity error'},'FontSize',14,'FontName','Arial');
box off
print -dpng error_subjects


% MAGNITUDE

subplot(2,2,1);
for s=1:round(length(subjects)/4)
plot(0:1:30,histc(sub_sacc_mag{s},0:1:30),'LineWidth',2,'Color',[randi([1 100],1,1)/100, randi([1 100],1,1)/100, randi([1 100],1,1)/100]); hold on
end
xlabel({'Saccade magnitude'},'FontSize',14,'FontName','Arial');
ylabel({'N'},'FontSize',14,'FontName','Arial');
box off;
 set(gca,'YLim',[0 300]);
subplot(2,2,2);
for s=round(round(length(subjects)/4)):round(length(subjects)/2)
plot(0:1:30,histc(sub_sacc_mag{s},0:1:30),'LineWidth',2,'Color',[randi([1 100],1,1)/100, randi([1 100],1,1)/100, randi([1 100],1,1)/100]); hold on
end
xlabel({'Saccade magnitude'},'FontSize',14,'FontName','Arial');
ylabel({'N'},'FontSize',14,'FontName','Arial');
box off;
 set(gca,'YLim',[0 300]);
subplot(2,2,3);
for s=round(length(subjects)/2):round(3*length(subjects)/4)
plot(0:1:30,histc(sub_sacc_mag{s},0:1:30),'LineWidth',2,'Color',[randi([1 100],1,1)/100, randi([1 100],1,1)/100, randi([1 100],1,1)/100]); hold on
end
xlabel({'Saccade magnitude'},'FontSize',14,'FontName','Arial');
ylabel({'N'},'FontSize',14,'FontName','Arial');
box off;
 set(gca,'YLim',[0 300]);
subplot(2,2,4);
for s=round(3*length(subjects)/4):length(subjects)
plot(0:1:30,histc(sub_sacc_mag{s},0:1:30),'LineWidth',2,'Color',[randi([1 100],1,1)/100, randi([1 100],1,1)/100, randi([1 100],1,1)/100]); hold on
end
xlabel({'Saccade magnitude'},'FontSize',14,'FontName','Arial');
ylabel({'N'},'FontSize',14,'FontName','Arial');
box off;
 set(gca,'YLim',[0 300]);

print -dpng magnitude_subjects



% Study based on cuts and using median (~=5 cuts)
% cuts_low = find(cuts < median(cuts)); %99
% cuts_high = find(cuts >median(cuts));% 85
disp('CUTS');
cuts_low_sacc_cvn = [];
cuts_high_sacc_cvn = [];
cuts_low_sacc_err = [];
cuts_high_sacc_err = [];
cuts_low_sacc_mag = [];
cuts_high_sacc_mag = [];
cuts_med_sacc_mag = [];
cuts_med_sacc_cvn = [];
cuts_med_sacc_err = [];
maglimit = 20;
for j=1:length(saccadeIndex)
    if cuts(saccadeIndex(j).videoNumber) < 4
        cuts_low_sacc_cvn = vertcat(cuts_low_sacc_cvn, saccadeIndex(j).curvature);
        cuts_low_sacc_err = vertcat(cuts_low_sacc_err, saccadeIndex(j).err);
        cuts_low_sacc_mag = vertcat(cuts_low_sacc_mag, saccadeIndex(j).length);
        
    else  if  cuts(saccadeIndex(j).videoNumber) >5
            cuts_high_sacc_cvn = vertcat(cuts_high_sacc_cvn, saccadeIndex(j).curvature);
            cuts_high_sacc_err = vertcat(cuts_high_sacc_err, saccadeIndex(j).err);
            cuts_high_sacc_mag = vertcat(cuts_high_sacc_mag, saccadeIndex(j).length);
            
        else
            cuts_med_sacc_cvn = vertcat(cuts_med_sacc_cvn, saccadeIndex(j).curvature);
            cuts_med_sacc_err = vertcat(cuts_med_sacc_err, saccadeIndex(j).err);
            cuts_med_sacc_mag = vertcat(cuts_med_sacc_mag, saccadeIndex(j).length);
        end
    end
end
datalow_cvn= histc(cuts_low_sacc_cvn,0:0.01:0.25)./ sum (histc(cuts_low_sacc_cvn,0:0.01:0.25));
datahigh_cvn= histc(cuts_high_sacc_cvn,0:0.01:0.25)./ sum (histc(cuts_high_sacc_cvn,0:0.01:0.25));
datamed_cvn= histc(cuts_med_sacc_cvn,0:0.01:0.25)./ sum (histc(cuts_med_sacc_cvn,0:0.01:0.25));
datalow_err = histc(cuts_low_sacc_err,0:0.0005:0.01)./ sum (histc(cuts_low_sacc_err,0:0.0005:0.01));
datahigh_err = histc(cuts_high_sacc_err,0:0.0005:0.01)./ sum (histc(cuts_high_sacc_err,0:0.0005:0.01));
datamed_err = histc(cuts_med_sacc_err,0:0.0005:0.01)./ sum (histc(cuts_med_sacc_err,0:0.0005:0.01));
datalow_mag= histc(cuts_low_sacc_mag,0:1:maglimit)./ sum (histc(cuts_low_sacc_mag,0:1:maglimit));
datahigh_mag= histc(cuts_high_sacc_mag,0:1:maglimit)./ sum (histc(cuts_high_sacc_mag,0:1:maglimit));
datamed_mag= histc(cuts_med_sacc_mag,0:1:maglimit)./ sum (histc(cuts_med_sacc_mag,0:1:maglimit));

[h_cvn p_cvn] = ttest(datalow_cvn, datahigh_cvn, 'alpha',0.05)
[h_err p_err] = ttest(datalow_err, datahigh_err, 'alpha',0.05)
figure1= figure;
subplot(2,1,1,'Parent',figure1,'FontSize',16,'FontName','Arial');
plot(0:0.01:0.25,datalow_cvn,'b','LineWidth', 2); hold on
plot(0:0.01:0.25,datahigh_cvn,'r','LineWidth', 2); hold on
plot(0:0.01:0.25,datamed_cvn,'g','LineWidth', 2); hold on
box off;
xlabel({'Saccade Curvature'},'FontSize',18,'FontName','Arial');
ylabel({'Frequency'},'FontSize',18,'FontName','Arial');
legend({'Low #Cuts','Medium #Cuts','High #Cuts'},'FontSize',16,'FontName','Arial');
%text(0.35,0.4,['H =' num2str(h_cvn) ' P-value =' num2str(p_cvn)]);
legend boxoff
subplot(2,1,2,'Parent',figure1,'FontSize',16,'FontName','Arial');
plot(0:0.0005:0.01,datalow_err,'b','LineWidth', 2); hold on;
plot(0:0.0005:0.01,datahigh_err,'r','LineWidth', 2); hold on;
plot(0:0.0005:0.01,datamed_err,'g','LineWidth', 2); hold on;
%text(0.007,0.4,['H =' num2str(h_err) ' P-value =' num2str(p_err)]);
box off;
ylabel({'Frequency'},'FontSize',18,'FontName','Arial');
xlabel({'Saccade Velocity error'},'FontSize',18,'FontName','Arial');
legend({'Low #Cuts','Medium #Cuts','High #Cuts'},'FontSize',16,'FontName','Arial');
legend boxoff
% subplot(,1,3,'Parent',figure1,'FontSize',16,'FontName','Arial');
% plot(0:1:maglimit,datalow_mag,'b','LineWidth', 2); hold on;
% plot(0:1:maglimit,datahigh_mag,'r','LineWidth', 2); hold on;
% plot(0:1:maglimit,datamed_mag,'g','LineWidth', 2); hold on;
% %text(0.007,0.4,['H =' num2str(h_err) ' P-value =' num2str(p_err)]);
% box off;
% ylabel({'N (Normalized)'},'FontSize',14,'FontName','Arial');
% xlabel({'Saccade Magnitude'},'FontSize',14,'FontName','Arial');
% legend({'Horizontal','Vertical','Oblique'},'FontSize',16,'FontName','Arial');
% legend boxoff
print -dpng cuts

%%%%
%old version with only 2 states for cuts

for j=1:length(saccadeIndex)
    if cuts(saccadeIndex(j).videoNumber) < median(cuts)
        cuts_low_sacc_cvn = vertcat(cuts_low_sacc_cvn, saccadeIndex(j).curvature);
        cuts_low_sacc_err = vertcat(cuts_low_sacc_err, saccadeIndex(j).err);
        cuts_low_sacc_mag = vertcat(cuts_low_sacc_mag, saccadeIndex(j).length);
    else if cuts(saccadeIndex(j).videoNumber) > median(cuts)
        cuts_high_sacc_cvn = vertcat(cuts_high_sacc_cvn, saccadeIndex(j).curvature);        
        cuts_high_sacc_err = vertcat(cuts_high_sacc_err, saccadeIndex(j).err);
         cuts_high_sacc_mag = vertcat(cuts_high_sacc_mag, saccadeIndex(j).length);
        end
    end
end

plot_histograms_odd_saccades(cuts_low_sacc_cvn, cuts_high_sacc_cvn, cuts_low_sacc_err, cuts_high_sacc_err, cuts_low_sacc_mag, cuts_high_sacc_mag, 0, 'Cuts');
% 
% datalow_cvn= histc(cuts_low_sacc_cvn,0:0.01:0.3)./ sum (histc(cuts_low_sacc_cvn,0:0.01:0.3));
% datahigh_cvn= histc(cuts_high_sacc_cvn,0:0.01:0.3)/sum (histc(cuts_low_sacc_cvn,0:0.01:0.3));
% datalow_err = histc(cuts_low_sacc_err,0:0.0005:0.01)/ sum (histc(cuts_low_sacc_err,0:0.0005:0.01));
% datahigh_err = histc(cuts_high_sacc_err,0:0.0005:0.01)/ sum (histc(cuts_high_sacc_err,0:0.0005:0.01));
% [h_cvn p_cvn] = ttest(datalow_cvn, datahigh_cvn, 'alpha',0.05);
% [h_err p_err] = ttest(datalow_err, datahigh_err, 'alpha',0.05);
% figure1= figure;
% subplot(2,1,1,'Parent',figure1,'FontSize',14,'FontName','Arial');
% plot(0:0.01:0.3,datalow_cvn,'b','LineWidth', 2); hold on
% plot(0:0.01:0.3,datahigh_cvn,'r','LineWidth', 2); hold on
% box off
% xlabel({'Saccade Curvature'},'FontSize',14,'FontName','Arial');
% ylabel({'N (Normalized)'},'FontSize',14,'FontName','Arial');
% legend({'Low Number of Cuts','High number of Cuts'},'FontSize',16,'FontName','Arial');
% legend boxoff
% % text(0.35,0.4,['H =' num2str(h_cvn) ' P-value =' num2str(p_cvn)]);
% subplot(2,1,2,'Parent',figure1,'FontSize',14,'FontName','Arial');
% plot(0:0.0005:0.01,datalow_err,'b','LineWidth', 2); hold on;
% plot(0:0.0005:0.01,datahigh_err,'r','LineWidth', 2); hold on;
% legend({'Low Number of Cuts','High number of Cuts'},'FontSize',16,'FontName','Arial');
% legend boxoff
% % text(0.007,0.4,['H =' num2str(h_err) ' P-value =' num2str(p_err)]);
% ylabel({'N (Normalized)'},'FontSize',14,'FontName','Arial');
% xlabel({'Saccade Velocity error'},'FontSize',14,'FontName','Arial');
% box off;
% print -dpng cuts

data1(:,1) = 1:length(datalow_cvn);
data1(:,2) = 1;
data1(:,3) = datalow_cvn;
data2(:,1) = 1:length(datalow_cvn);
data2(:,2) = 2;
data2(:,3) = datahigh_cvn;
cuts_cvn = [ data1; data2];
dlmwrite('cuts_cvn.txt',cuts_cvn,'delimiter',',')


%% Study based on order
disp('Order');
order_low_sacc_cvn = []; %88
order_high_sacc_cvn = []; % 118
order_low_sacc_err = [];
order_high_sacc_err = [];
order_low_sacc_mag = [];
order_high_sacc_mag = [];


for j=1:length(saccadeIndex)
    if saccadeIndex(j).order < 24
        order_low_sacc_cvn = vertcat(order_low_sacc_cvn, saccadeIndex(j).curvature);
        order_low_sacc_err = vertcat(order_low_sacc_err, saccadeIndex(j).err);
        order_low_sacc_mag = vertcat(order_low_sacc_mag, saccadeIndex(j).length);
   
    else
        order_high_sacc_cvn = vertcat(order_high_sacc_cvn, saccadeIndex(j).curvature);        
        order_high_sacc_err = vertcat(order_high_sacc_err, saccadeIndex(j).err);
        order_high_sacc_mag = vertcat(order_high_sacc_mag, saccadeIndex(j).length);
   
    end
end

mean(order_low_sacc_cvn);
mean(order_high_sacc_cvn);
plot_histograms_odd_saccades(order_low_sacc_cvn, order_high_sacc_cvn, order_low_sacc_err, order_high_sacc_err, order_low_sacc_mag, order_high_sacc_mag, 0, 'order');



%% Study based on nature
disp('Nature');
nature_low_sacc_cvn = []; %88
nature_high_sacc_cvn = []; % 118
nature_low_sacc_err = [];
nature_high_sacc_err = [];
nature_low_sacc_mag = [];
nature_high_sacc_mag = [];

for j=1:length(saccadeIndex)
    if nature05(saccadeIndex(j).videoNumber) == 0
        nature_low_sacc_cvn = vertcat(nature_low_sacc_cvn, saccadeIndex(j).curvature);
        nature_low_sacc_err = vertcat(nature_low_sacc_err, saccadeIndex(j).err);
          nature_low_sacc_mag = vertcat(nature_low_sacc_mag, saccadeIndex(j).length);
   
    else
        nature_high_sacc_cvn = vertcat(nature_high_sacc_cvn, saccadeIndex(j).curvature);        
        nature_high_sacc_err = vertcat(nature_high_sacc_err, saccadeIndex(j).err);
        nature_high_sacc_mag = vertcat(nature_high_sacc_mag, saccadeIndex(j).length);
   
    end
end

plot_histograms_odd_saccades(nature_low_sacc_cvn, nature_high_sacc_cvn, nature_low_sacc_err, nature_high_sacc_err, nature_low_sacc_mag, nature_high_sacc_mag, 0, 'Nature');



%% Study based on light
% discarding condition 3
light_low_sacc_cvn = []; %67
light_high_sacc_cvn = []; % 51
light_low_sacc_err = [];
light_high_sacc_err = [];
light_low_sacc_mag = [];
light_high_sacc_mag = [];
disp('Light');
for j=1:length(saccadeIndex)
    if light05(saccadeIndex(j).videoNumber) < median(light05)
        light_low_sacc_cvn = vertcat(light_low_sacc_cvn, saccadeIndex(j).curvature);
        light_low_sacc_err = vertcat(light_low_sacc_err, saccadeIndex(j).err);
         light_low_sacc_mag = vertcat(light_low_sacc_mag, saccadeIndex(j).length);
    else  if light05(saccadeIndex(j).videoNumber) > median(light05)
        light_high_sacc_cvn = vertcat(light_high_sacc_cvn, saccadeIndex(j).curvature);        
        light_high_sacc_err = vertcat(light_high_sacc_err, saccadeIndex(j).err);
           light_high_sacc_mag = vertcat(light_high_sacc_mag, saccadeIndex(j).length);
        end
    end
end

plot_histograms_odd_saccades(light_low_sacc_cvn, light_high_sacc_cvn, light_low_sacc_err, light_high_sacc_err,light_low_sacc_mag, light_high_sacc_mag, 0, 'Light');



%% Study based on AudioInfo
% discarding condition 3
audio_low_sacc_cvn = []; %83
audio_high_sacc_cvn = []; % 63
audio_low_sacc_err = [];
audio_high_sacc_err = [];
audio_low_sacc_mag = [];
audio_high_sacc_mag = [];
disp('Audio');
for j=1:length(saccadeIndex)
    if audInfo(saccadeIndex(j).videoNumber) < nanmedian(audInfo)
        audio_low_sacc_cvn = vertcat(audio_low_sacc_cvn, saccadeIndex(j).curvature);
        audio_low_sacc_err = vertcat(audio_low_sacc_err, saccadeIndex(j).err);
         audio_low_sacc_mag = vertcat(audio_low_sacc_mag, saccadeIndex(j).length);
    else  if audInfo(saccadeIndex(j).videoNumber) > nanmedian(audInfo)
        audio_high_sacc_cvn = vertcat(audio_high_sacc_cvn, saccadeIndex(j).curvature);        
        audio_high_sacc_err = vertcat(audio_high_sacc_err, saccadeIndex(j).err);
        audio_high_sacc_mag = vertcat(audio_high_sacc_mag, saccadeIndex(j).length);
        end
    end
end

plot_histograms_odd_saccades(audio_low_sacc_cvn, audio_high_sacc_cvn, audio_low_sacc_err, audio_high_sacc_err, audio_low_sacc_mag, audio_high_sacc_mag,0, 'AudioInfo');



%% Study based on Faces
% discarding 2 faces
faces_low_sacc_cvn = []; %82
faces_high_sacc_cvn = []; % 71
faces_low_sacc_err = [];
faces_high_sacc_err = [];
faces_low_sacc_mag = [];
faces_high_sacc_mag = [];
disp('Faces');
for j=1:length(saccadeIndex)
    if faces05(saccadeIndex(j).videoNumber) < nanmedian(faces05)
        faces_low_sacc_cvn = vertcat(faces_low_sacc_cvn, saccadeIndex(j).curvature);
        faces_low_sacc_err = vertcat(faces_low_sacc_err, saccadeIndex(j).err);
         faces_low_sacc_mag = vertcat(faces_low_sacc_mag, saccadeIndex(j).length);
    else  if faces05(saccadeIndex(j).videoNumber) > nanmedian(faces05)
        faces_high_sacc_cvn = vertcat(faces_high_sacc_cvn, saccadeIndex(j).curvature);        
        faces_high_sacc_err = vertcat(faces_high_sacc_err, saccadeIndex(j).err);
         faces_high_sacc_mag = vertcat(faces_high_sacc_mag, saccadeIndex(j).length);
        end
    end
end

plot_histograms_odd_saccades(faces_low_sacc_cvn, faces_high_sacc_cvn, faces_low_sacc_err, faces_high_sacc_err, faces_low_sacc_mag, faces_high_sacc_mag,0, 'Face');



%% Environment Indoor/outdoor

envi_low_sacc_cvn = []; % indoor %92
envi_high_sacc_cvn = []; % outdoor %98
envi_low_sacc_err = [];
envi_high_sacc_err = [];
envi_low_sacc_mag = [];
envi_high_sacc_mag = [];
disp('envi');
for j=1:length(saccadeIndex)
    if strcmp(environment(saccadeIndex(j).videoNumber),'indoor')
        envi_low_sacc_cvn = vertcat(envi_low_sacc_cvn, saccadeIndex(j).curvature);
        envi_low_sacc_err = vertcat(envi_low_sacc_err, saccadeIndex(j).err);
          envi_low_sacc_mag = vertcat(envi_low_sacc_mag, saccadeIndex(j).length);
    else  
        envi_high_sacc_cvn = vertcat(envi_high_sacc_cvn, saccadeIndex(j).curvature);        
        envi_high_sacc_err = vertcat(envi_high_sacc_err, saccadeIndex(j).err);
         envi_high_sacc_mag = vertcat(envi_high_sacc_mag, saccadeIndex(j).length);
        
    end
end

plot_histograms_odd_saccades(envi_low_sacc_cvn, envi_high_sacc_cvn, envi_low_sacc_err, envi_high_sacc_err, envi_low_sacc_mag, envi_high_sacc_mag, 0, 'Environment');



%% According to direction
dire_low_sacc_cvn = []; % indoor %92
dire_high_sacc_cvn = []; % outdoor %98
dire_low_sacc_err = [];
dire_high_sacc_err = [];
dire_med_sacc_cvn = [];
dire_med_sacc_err = [];
dire_med_sacc_mag = [];
dire_high_sacc_mag = [];
dire_low_sacc_mag = [];
maglimit = 20;
disp('dire');
% This set is for 5 degrees in horizontal
for j=1:length(saccadeIndex)
    if (cos(deg2rad(saccadeIndex(j).direction)) > 0.99) || (cos(deg2rad(saccadeIndex(j).direction)) < -0.99)
        dire_low_sacc_cvn = vertcat(dire_low_sacc_cvn, saccadeIndex(j).curvature);
        dire_low_sacc_err = vertcat(dire_low_sacc_err, saccadeIndex(j).err);
        dire_low_sacc_mag = vertcat(dire_low_sacc_mag, saccadeIndex(j).length);
       
    else  if  (cos(deg2rad(saccadeIndex(j).direction)) > - 0.3827) && (cos(deg2rad(saccadeIndex(j).direction)) < 0.3827)
        dire_high_sacc_cvn = vertcat(dire_high_sacc_cvn, saccadeIndex(j).curvature);        
        dire_high_sacc_err = vertcat(dire_high_sacc_err, saccadeIndex(j).err);
        dire_high_sacc_mag = vertcat(dire_high_sacc_mag, saccadeIndex(j).length);
        
        else
            dire_med_sacc_cvn = vertcat(dire_med_sacc_cvn, saccadeIndex(j).curvature);        
            dire_med_sacc_err = vertcat(dire_med_sacc_err, saccadeIndex(j).err); 
            dire_med_sacc_mag = vertcat(dire_med_sacc_mag, saccadeIndex(j).length);
        end
    end
end
% for j=1:length(saccadeIndex)
%     if (cos(deg2rad(saccadeIndex(j).direction)) > 0.92) || (cos(deg2rad(saccadeIndex(j).direction)) < -0.92)
%         dire_low_sacc_cvn = vertcat(dire_low_sacc_cvn, saccadeIndex(j).curvature);
%         dire_low_sacc_err = vertcat(dire_low_sacc_err, saccadeIndex(j).err);
%         dire_low_sacc_mag = vertcat(dire_low_sacc_mag, saccadeIndex(j).length);
%        
%     else  if  (cos(deg2rad(saccadeIndex(j).direction)) > - 0.3827) && (cos(deg2rad(saccadeIndex(j).direction)) < 0.3827)
%         dire_high_sacc_cvn = vertcat(dire_high_sacc_cvn, saccadeIndex(j).curvature);        
%         dire_high_sacc_err = vertcat(dire_high_sacc_err, saccadeIndex(j).err);
%         dire_high_sacc_mag = vertcat(dire_high_sacc_mag, saccadeIndex(j).length);
%         
%         else
%             dire_med_sacc_cvn = vertcat(dire_med_sacc_cvn, saccadeIndex(j).curvature);        
%             dire_med_sacc_err = vertcat(dire_med_sacc_err, saccadeIndex(j).err); 
%             dire_med_sacc_mag = vertcat(dire_med_sacc_mag, saccadeIndex(j).length);
%         end
%     end
% end

datalow_cvn= histc(dire_low_sacc_cvn,0:0.01:0.25)./ sum (histc(dire_low_sacc_cvn,0:0.01:0.25));
datahigh_cvn= histc(dire_high_sacc_cvn,0:0.01:0.25)./ sum (histc(dire_high_sacc_cvn,0:0.01:0.25));
datamed_cvn= histc(dire_med_sacc_cvn,0:0.01:0.25)./ sum (histc(dire_med_sacc_cvn,0:0.01:0.25));
datalow_err = histc(dire_low_sacc_err,0:0.0005:0.01)./ sum (histc(dire_low_sacc_err,0:0.0005:0.01));
datahigh_err = histc(dire_high_sacc_err,0:0.0005:0.01)./ sum (histc(dire_high_sacc_err,0:0.0005:0.01));
datamed_err = histc(dire_med_sacc_err,0:0.0005:0.01)./ sum (histc(dire_med_sacc_err,0:0.0005:0.01));
datalow_mag= histc(dire_low_sacc_mag,0:1:maglimit)./ sum (histc(dire_low_sacc_mag,0:1:maglimit));
datahigh_mag= histc(dire_high_sacc_mag,0:1:maglimit)./ sum (histc(dire_high_sacc_mag,0:1:maglimit));
datamed_mag= histc(dire_med_sacc_mag,0:1:maglimit)./ sum (histc(dire_med_sacc_mag,0:1:maglimit));

[h_cvn p_cvn] = ttest(datalow_cvn, datahigh_cvn, 'alpha',0.05)
[h_err p_err] = ttest(datalow_err, datahigh_err, 'alpha',0.05)
figure1= figure;
subplot(3,1,1,'Parent',figure1,'FontSize',16,'FontName','Arial');
plot(0:0.01:0.25,datalow_cvn,'b','LineWidth', 2); hold on
plot(0:0.01:0.25,datahigh_cvn,'r','LineWidth', 2); hold on
plot(0:0.01:0.25,datamed_cvn,'g','LineWidth', 2); hold on
box off;
xlabel({'Saccade Curvature'},'FontSize',18,'FontName','Arial');
ylabel({'N (Normalized)'},'FontSize',18,'FontName','Arial');
legend({'Horizontal','Vertical','Oblique'},'FontSize',16,'FontName','Arial');
%text(0.35,0.4,['H =' num2str(h_cvn) ' P-value =' num2str(p_cvn)]);
legend boxoff
subplot(3,1,2,'Parent',figure1,'FontSize',16,'FontName','Arial');
plot(0:0.0005:0.01,datalow_err,'b','LineWidth', 2); hold on;
plot(0:0.0005:0.01,datahigh_err,'r','LineWidth', 2); hold on;
plot(0:0.0005:0.01,datamed_err,'g','LineWidth', 2); hold on;
%text(0.007,0.4,['H =' num2str(h_err) ' P-value =' num2str(p_err)]);
box off;
ylabel({'N (Normalized)'},'FontSize',18,'FontName','Arial');
xlabel({'Saccade Velocity error'},'FontSize',18,'FontName','Arial');
legend({'Horizontal','Vertical','Oblique'},'FontSize',16,'FontName','Arial');
legend boxoff
subplot(3,1,3,'Parent',figure1,'FontSize',16,'FontName','Arial');
plot(0:1:maglimit,datalow_mag,'b','LineWidth', 2); hold on;
plot(0:1:maglimit,datahigh_mag,'r','LineWidth', 2); hold on;
plot(0:1:maglimit,datamed_mag,'g','LineWidth', 2); hold on;
%text(0.007,0.4,['H =' num2str(h_err) ' P-value =' num2str(p_err)]);
box off;
ylabel({'N (Normalized)'},'FontSize',14,'FontName','Arial');
xlabel({'Saccade Magnitude'},'FontSize',14,'FontName','Arial');
legend({'Horizontal','Vertical','Oblique'},'FontSize',16,'FontName','Arial');
legend boxoff
print -dpng direction

%% According to age
age_low_sacc_cvn = []; % indoor %92
age_high_sacc_cvn = []; % outdoor %98
age_low_sacc_err = [];
age_high_sacc_err = [];
age_med_sacc_cvn = [];
age_med_sacc_err = [];
age_low_sacc_mag = [];
age_med_sacc_mag = [];
age_high_sacc_mag = [];
disp('age');
for j=1:length(saccadeIndex)
    if (saccadeIndex(j).age < y(1))
        age_low_sacc_cvn = vertcat(age_low_sacc_cvn, saccadeIndex(j).curvature);
        age_low_sacc_err = vertcat(age_low_sacc_err, saccadeIndex(j).err);
        age_low_sacc_mag = vertcat(age_low_sacc_mag, saccadeIndex(j).length);
    else  if  (saccadeIndex(j).age > y(2))
        age_high_sacc_cvn = vertcat(age_high_sacc_cvn, saccadeIndex(j).curvature);        
        age_high_sacc_err = vertcat(age_high_sacc_err, saccadeIndex(j).err);
        age_high_sacc_mag = vertcat(age_high_sacc_mag, saccadeIndex(j).length);
        
        else
            age_med_sacc_cvn = vertcat(age_med_sacc_cvn, saccadeIndex(j).curvature);        
            age_med_sacc_err = vertcat(age_med_sacc_err, saccadeIndex(j).err); 
            age_med_sacc_mag = vertcat(age_med_sacc_mag, saccadeIndex(j).length);
        end
    end
end
datalow_cvn= histc(age_low_sacc_cvn,0:0.01:0.25)./ sum (histc(age_low_sacc_cvn,0:0.01:0.25));
datahigh_cvn= histc(age_high_sacc_cvn,0:0.01:0.25)./ sum (histc(age_high_sacc_cvn,0:0.01:0.25));
datamed_cvn= histc(age_med_sacc_cvn,0:0.01:0.25)./ sum (histc(age_med_sacc_cvn,0:0.01:0.25));
datalow_err = histc(age_low_sacc_err,0:0.0005:0.01)./ sum (histc(age_low_sacc_err,0:0.0005:0.01));
datahigh_err = histc(age_high_sacc_err,0:0.0005:0.01)./ sum (histc(age_high_sacc_err,0:0.0005:0.01));
datamed_err = histc(age_med_sacc_err,0:0.0005:0.01)./ sum (histc(age_med_sacc_err,0:0.0005:0.01));
datalow_mag = histc(age_low_sacc_mag,0:1:40)./ sum (histc(age_low_sacc_mag,0:1:40));
datahigh_mag = histc(age_high_sacc_mag,0:1:40)./ sum (histc(age_high_sacc_mag,0:1:40));
datamed_mag = histc(age_med_sacc_mag,0:1:40)./ sum (histc(age_med_sacc_mag,0:1:40));
% [h_cvn p_cvn] = ttest(datalow_cvn, datahigh_cvn, 'alpha',0.05)
% [h_err p_err] = ttest(datalow_err, datahigh_err, 'alpha',0.05)
figure1= figure;
subplot(3,1,1,'Parent',figure1,'FontSize',14,'FontName','Arial');
plot(0:0.01:0.25,datalow_cvn,'b','LineWidth', 2); hold on
plot(0:0.01:0.25,datahigh_cvn,'r','LineWidth', 2); hold on
plot(0:0.01:0.25,datamed_cvn,'g','LineWidth', 2); hold on
box off;
xlabel({'Saccade Curvature'},'FontSize',14,'FontName','Arial');
ylabel({'N (Normalized)'},'FontSize',14,'FontName','Arial');
legend({'<43','>61','43-61'},'FontSize',16,'FontName','Arial');
%text(0.35,0.4,['H =' num2str(h_cvn) ' P-value =' num2str(p_cvn)]);
legend boxoff
subplot(3,1,2,'Parent',figure1,'FontSize',14,'FontName','Arial');
plot(0:0.0005:0.01,datalow_err,'b','LineWidth', 2); hold on;
plot(0:0.0005:0.01,datahigh_err,'r','LineWidth', 2); hold on;
plot(0:0.0005:0.01,datamed_err,'g','LineWidth', 2); hold on;
%text(0.007,0.4,['H =' num2str(h_err) ' P-value =' num2str(p_err)]);
box off;
ylabel({'N (Normalized)'},'FontSize',14,'FontName','Arial');
xlabel({'Saccade Velocity error'},'FontSize',14,'FontName','Arial');
legend({'<43','>61','43-61'},'FontSize',16,'FontName','Arial');
legend boxoff
subplot(3,1,3,'Parent',figure1,'FontSize',14,'FontName','Arial');
plot(0:1:40,datalow_mag,'b','LineWidth', 2); hold on
plot(0:1:40,datahigh_mag,'r','LineWidth', 2); hold on
plot(0:1:40,datamed_mag,'g','LineWidth', 2); hold on
box off;
xlabel({'Saccade Magnitude'},'FontSize',14,'FontName','Arial');
ylabel({'N (Normalized)'},'FontSize',14,'FontName','Arial');
legend({'<43','>61','43-61'},'FontSize',16,'FontName','Arial');
%text(0.35,0.4,['H =' num2str(h_cvn) ' P-value =' num2str(p_cvn)]);
legend boxoff
print -dpng age

%% Plots of Curvatures and error fits

    cvrange=[.09 .16]%[.17 .3];
    cvind=find(medcvn>=cvrange(1) & medcvn<cvrange(2))
    figure;set(gcf,'color',[1 1 1])
    for i=1:length(cvind)
        x=sac{cvind(i)}(:,2)-sac{cvind(i)}(1,2);
        y=sac{cvind(i)}(:,1)-sac{cvind(i)}(1,1);
        plot(x,y,'o-');
        set(gca,'xlim',[-40 40],'ylim',[-40 40],...
            'xtick',[-40:5:40],'ytick',[-40:5:40]);
        xlabel('horz (deg)');
        ylabel('vert (deg)');
        text(-35,35,['curvature =' num2str(medcvn(cvind(i))) ' saccade # ' num2str(cvind(i)) '  Error:' num2str(err(cvind(i))) ], 'FontSize', 18);
        %     set(gca,'xlim',[-500 500],'ylim',[-500 500],...
        %             'xtick',[-500:100:500],'ytick',[-500:100:500])
        %     xlabel('horz (pix)')
        %     ylabel('vert (pix)')
        text(-450,450,['curvature =' num2str(medcvn(cvind(i))) ' saccade # ' num2str(cvind(i))  ]);
        axis square
        grid on
        pause
    end
    
    %% plot saccade and velocity fit
    %load allSaccades3
    %load allSaccades3_velocity
   % sind=find(medcvn<.45 & err<.004); % It was 0.045 before
    %sind=[5221 7593 7594 10073 11302 11670 15701 15752 16443]; % saccades that change direction
   % sind=find(medcvn>.20 & err>.004); % It was 0.045 before
     sind=find(err>.005); % It was 0.045 before
    
    figure;set(gcf,'color',[1 1 1])                            % indices from allSaccades2
    for i=1:length(sind)
        %plot saccade
        subplot(1,2,1)
        x=sac{sind(i)}(:,2)-sac{sind(i)}(1,2);
        y=sac{sind(i)}(:,1)-sac{sind(i)}(1,1);
        plot(x,y,'o-')
        %     set(gca,'xlim',[-700 700],'ylim',[-500 500],...
        %             'xtick',[-700:100:700],'ytick',[-500:100:500])
        set(gca,'xlim',[-20 20],'ylim',[-20 20],...
            'xtick',[-20:2:20],'ytick',[-20:2:20])
        axis square
        grid on
        xlabel('horz (deg)')
        ylabel('vert (deg)')
        text(-18,18,['curvature =' num2str(medcvn(sind(i)))])
        
        % plot velocity and fit
        subplot(1,2,2)
        plot(vel{sind(i)},'bo-','markersize',12);hold on
        plot(velfit{sind(i)},'r','linewidth',2);hold off
        set(gca,'ylim',[0 1],'ytick',0:.1:1)
        xlabel('t (ms)')
        ylabel('speed (deg/ms)')
        title(['err = ' num2str(err(sind(i))) '   saccade # ' num2str(sind(i))])
        pause
    end
    
    avgcur= [];
    avgerr = [];

for i=1:60
    ind = find(order==i);
    avgcur(i) = median(medcvn(ind));
    avgerr(i) = median(err(ind));
end

figure;plot(1:1:60,avgcur,'r','Marker','o', 'LineWidth', 2);
set(gca,'FontSize',16);
xlabel('Time-on-task (30-seconds bins)','FontSize', 18, 'FontName', 'Arial');
ylabel('Curvature','FontSize', 18, 'FontName', 'Arial');
box off
hold on
 z = smooth(avgcur);
plot(1:1:60,z,'k', 'LineWidth', 4);

