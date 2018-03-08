%This scripts reads from the re-constructed undersample data
clear all


% set path to eye tracker data, load video numbers and subject indices
etPath='/Users/FranciscoCostela/Desktop/magnification/Videoclipeyetrackingdata';
etPathN='/Users/FranciscoCostela/Desktop/GazeContigentDisplay/matlab_undersample';

outPath = '/Users/FranciscoCostela/Desktop/GazeContigentDisplay/SaccadeDB/';

%'/Users/FranciscoCostela/Desktop/AckermannArchive/Video clip eyetracking data';
load 'video number lookup.mat'
load baselineSubjectIndex

%etPath='/Users/johnackermann/Desktop/SchepensArchive/OfficeComp/Documents/Video clip eyetracking data';
%load '/Users/johnackermann/Desktop/SchepensArchive/OfficeComp/Documents/video number lookup.mat';
%load '/Users/johnackermann/Desktop/SchepensArchive/OfficeComp/Documents/baselineSubjectIndex';

% name for the output .mat file and struct index
filename='allSaccadesTV';
structname='saccadeIndexTV_plus16before';

% saccade struct index
k=0;
sac={};
sacNoFilter={};
tim={};
vnum=[];
trialnum=[];
subnum=[];

screenDist=100;        % distance in cm
screenSize=[60 34];   % [horz vert] screen dimensions in cm
screenRes=[2560 1440];% screen dims in pixs

% loop through each clip
for vidNum=1:max(videoNumbers);
    
    disp(vidNum);
    % get all trials that used this clip
    trials4Vid=find(videoNumbers==vidNum);
    
    % get the subject number for each trial
    subs4Vid=subjectIndex(trials4Vid); % see subjectList for sub IDs
    
    % get all eye data for this clip
    for i = 1:length(trials4Vid)
        if exist([etPathN filesep eyetrackFiles{trials4Vid(i)}])>0
            load([etPathN filesep eyetrackFiles{trials4Vid(i)}]);
        else
            load([etPath filesep eyetrackFiles{trials4Vid(i)}]);
            
        end
                        
       % load([etPath filesep eyetrackFiles{trials4Vid(i)}]);
        temp.x = eyetrackRecord.x;
        temp.y = eyetrackRecord.y;
        temp.t = eyetrackRecord.t;
        temp.missing = eyetrackRecord.missing;                                                                                                                                                                                                                                                     
        etData(i) = temp;
    end
    
    % loop through each trial
    edint=cell(length(etData),1);
    edintNoFilter=cell(length(etData),1);
    sacind=cell(length(etData),1);
    timint=cell(length(etData),1);
    trialind=zeros(length(etData),1);
    subind=zeros(length(etData),1);
    for i=1:length(etData)
        if length(etData(i).t)>19000 % make sure there's enough data
            disp(['ok with' eyetrackFiles{trials4Vid(i)}]);
            ed=[etData(i).x' etData(i).y'];
            time=etData(i).t;
            % convert to degrees
            eddeg=pix2deg(ed,screenDist,screenSize,screenRes);
            
            % index saccades
          %  [edint{i},sacind{i},timint{i},fixpos, edintNoFilter{i}]=offlineSaccadeDetect_nofilter(eddeg,time); % allSaccades3 uses offlineSaccadeDetect2 _nofilter
            [edint{i},sacind{i},timint{i},fixpos, edintNoFilter{i}]=offlineSaccadeDetect_noSmooth(eddeg,time); % allSaccades3 uses offlineSaccadeDetect2 _nofilter
            
          % index trial
            trialind(i,1)=trials4Vid(i);
            % index subject
            subind(i,1)=subs4Vid(i);
        else
            disp(['not enough data with' eyetrackFiles{trials4Vid(i)}]); 
        end
    end
    
    % extract the saccades
    for sub=1:length(edint)
        insac=false;
        for i=1:length(sacind{sub})
            if sacind{sub}(i)==1 && ~insac
                k=k+1;
                sac{k}=edint{sub}(i,:); % saccade coordinates
                sacNoFilter{k}=edintNoFilter{sub}(max(1,i-16):i,:); % saccade coordinates
                tim{k}=timint{sub}(i);  % get clock time of the saccade
                vnum(k)=vidNum;         % the clip number 
                trialnum(k)=trialind(sub);  % the trial it's derived from (1:2681)
                subnum(k)=subind(sub);
                insac=true;
                continue
            elseif sacind{sub}(i)==1 && insac
                sac{k}(end+1,:)=edint{sub}(i,:);
                sacNoFilter{k}(end+1,:)=edintNoFilter{sub}(i,:);
                tim{k}(end+1,1)=timint{sub}(i);
                continue
            elseif sacind{sub}(i)==0 && insac
                insac=false;
                continue
            elseif sacind{sub}(i)==0 && ~insac
                continue
            end
        end
    end 
    eyeDataInterp{vidNum}=edint; % keep the full, interpolated data set
    trials{vidNum}=trials4Vid;   % corresponding trials
    subjects{vidNum}=subs4Vid;   % and subject indices
    clear etData
end

% get curvature
medcvn=zeros(length(sac),1);
for i=1:length(sac)
    X=sac{i}(:,1);
    Y=sac{i}(:,2); 
    [cv,cvn]=curvature(X,Y);
    medcvn(i)=median(cvn);
end

% get length
saclength=zeros(length(sac),1);
for i=1:length(sac)
    X=sac{i}(:,1);
    Y=sac{i}(:,2);
    saclength(i)=norm([X(end) Y(end)]-[X(1) Y(1)]);
end

% get direction
direction=zeros(length(sac),1);
for sacindex=1:length(sac)
    XY=sac{sacindex};
    endpoint=XY(end,:)-XY(1,:);
    ang=acos(endpoint(1)/norm(endpoint));
    if endpoint(2)<0
       ang=2*pi-ang;
    end
    direction(sacindex)=360 - ang*180/pi;
end

% get velocity of each saccade and fit exponential velocity fn
vel=cell(length(sac),1);
velfit=cell(length(sac),1);
params=cell(length(sac),1);
err=zeros(length(sac),1);   
for sacindex=1:length(sac)
     if ~mod(sacindex,10000)
        disp(sacindex)
    end
    XY=sac{sacindex};
    vel{sacindex}=sqrt(sum(diff(XY).^2,2));
    [velfit{sacindex},params{sacindex},err(sacindex)]=fitVelocity(vel{sacindex},4);
end

% % save 
% %save(filename,'sac','tim','vnum','trialnum',...
%               'subnum','eyeDataInterp','trials','subjects',...
%               'medcvn','saclength','direction')
%save([filename '_velocity'],'vel','velfit','params','err')

% creat struct index
j = 1;
k=1;
saccadeIndex=struct;

for i=1:length(sac)
    
    % We apply some thresholds to get rid of weird saccades
    % 1. length range: from 0.5 to 40 deg (saccades lower than 0.5 are considered micro-saccades, larger than 40 probably unrealistic given the data collection setup)
    % 2. duration < 100 ms
    % 3. terminal velocity < 0.3 deg/ms (for removing chopped saccades, this is set about 10 times the velocity for detection of saccade termination in original algorithm)
    % 4. DISCARDED initial velocity < 0.075 deg/ms (for removing chopped saccades, this is set about 2.5 times the velocity for detection of saccade initiation in original algorithm)
    % 5. velocity at 25% of duration > 0.15*peak velocity (this is to remove the slow starting saccades with a "knee", as such saccades tend to have a uniform but low velocity during their initial phase that covers first quartile of duration).    
    cx = abs( diff(sacNoFilter{i}(:,1)));
    cy = abs( diff(sacNoFilter{i}(:,2)));
    
    if ~any(cx>4) && ~any(cy>4) && length(sac{i})>15 && saclength(i) <40 && saclength(i)>1 && length(sac{i})<150 && vel{i}(end)<0.3 &&  vel{i}(round(length(sac{i})/4)) > 0.15*max(vel{i})
%         && vel{i}(1)<0.075 
        saccadeIndex(j).saccade=sac{i};
        saccadeIndex(j).duration = length(sac{i});
        saccadeIndex(j).saccadeNoFilter=sacNoFilter{i};        
        saccadeIndex(j).videoNumber=vnum(i);
        saccadeIndex(j).trialNumber=trialnum(i);
        saccadeIndex(j).subNumber=subnum(i);
        saccadeIndex(j).elapsedtime=tim{i};
        saccadeIndex(j).curvature=medcvn(i);
        saccadeIndex(j).length=saclength(i);
        saccadeIndex(j).direction=direction(i);
        saccadeIndex(j).velocity.velocity=vel{i};
        saccadeIndex(j).velocity.velocityFit=velfit{i};
        saccadeIndex(j).velocity.params=params{i};
        saccadeIndex(j).velocity.err=err(i);
        saccadeIndex(j).pkvelocity = max(vel{i});
      %  saccadeIndex(j).pushed = pushnum(i);
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
         saccadeIndex(j).type = 1;
        j = j+1;
    end
    
    
    
    saccadeIndexAll(k).saccade=sac{i}; 
    saccadeIndexAll(k).saccadeNoFilter=sacNoFilter{i};        
    saccadeIndexAll(k).trialNumber=trialnum(i);
    saccadeIndexAll(k).subNumber=subnum(i);
    saccadeIndexAll(k).elapsedtime=tim{i};
    saccadeIndexAll(k).curvature=medcvn(i);
    saccadeIndexAll(k).length=saclength(i);
    saccadeIndexAll(k).direction=direction(i);
    saccadeIndexAll(k).velocity.velocity=vel{i};
    saccadeIndexAll(k).velocity.velocityFit=velfit{i};
    saccadeIndexAll(k).velocity.params=params{i};
    saccadeIndexAll(k).velocity.err=err(i);
    saccadeIndexAll(k).pkvelocity = max(vel{i});  
    % +/- 5 degrees for horizontal . Reference:
    %  Duke Elder, S. & Wybar, K. (1973). Ocular Motility and Stabismus. (3rd ed.). System of Ophthalmology  (Vol. VI). London: Henry Kimpton.
    if (cos(deg2rad(saccadeIndexAll(k).direction)) > 0.99) || (cos(deg2rad(saccadeIndexAll(k).direction)) < -0.99)
        saccadeIndexAll(k).horver = 1; % horizontal
    else  if  (cos(deg2rad(saccadeIndexAll(k).direction)) > - 0.3827) && (cos(deg2rad(saccadeIndexAll(k).direction)) < 0.3827)
            saccadeIndexAll(k).horver = 3; % vertical
        else
            saccadeIndexAll(k).horver = 2; % oblique
        end
    end
    k = k +1;
    
end

% 
% nu = 1;
% close all
% figure
% plot (etData(nu).x)
% hold on
% plot(etData(nu).y)
% hold on
% %xlim([0 6000]);
% plot(sacind{1}*2000,'r')
% 
% 
% for i=1:length(saccadeIndex)
%     plot(saccadeIndex(i).elapsedtime,etData(nu).x(saccadeIndex(i).elapsedtime), 'g','LineWidth',5); 
%     hold on
%     plot(saccadeIndex(i).elapsedtime,etData(nu).y(saccadeIndex(i).elapsedtime), 'g','LineWidth',5); 
% end
% xlim([0 1700])
% 
% figure
% for i=1:6
%     i
%     
%     subplot(3,2,i);
%   
%    % text(20,20,num2str(saccadeIndex(i).direction));
%     plot(saccadeIndex(i).saccade(:,1), saccadeIndex(i).saccade(:,2),'k','LineWidth',2);
%     hold on
%     plot(saccadeIndex(i).saccadeNoFilter(:,1), saccadeIndex(i).saccadeNoFilter(:,2),'b','LineWidth',2);
%     xlabel('Horizontal (deg)');
%     ylabel('Vertical (deg)');
%     title(['Saccade ' num2str(i)],'FontSize', 18); 
%       axis equal
% end
% 
% for i=1:length(saccadeIndexAll)
%     velini(i) = saccadeIndexAll(i).velocity.velocity(1);
%     magi(i) = saccadeIndexAll(i).length;
%     pkvel(i) = saccadeIndexAll(i).pkvelocity;
% end
% 
% a = find(pkvel<0.1);
% b = velini(a);
% for i=100:200
% i
% scatter(magi(i), pkvel(i),10);
% disp(velini(i));
% hold on
% pause;
% end

% % 
% 
% l40 = 0;
% l05 = 0;
% l100 = 0;
% v03 = 0;
% v0075 = 0;
% v15 = 0;
% 
% for i=1:length(sac)
%         
%     if saclength(i) >40
%         
%         l40 = l40 +1;
%     end
%     if saclength(i)<0.5
%         l05 = l05 +1;
%     end
%     if length(sac{i})>100
%         l100 = l100 +1;
%     end
%     if vel{i}(end)>0.3
%         v03 = v03 +1;
%     end
%     if vel{i}(1)>0.075
%         v0075 = v0075 +1;
%     end
%     if vel{i}(round(length(sac{i})/4)) < 0.15*max(vel{i})
%         v15 = v15 +1;
%     end
%     
%     
% end
% 
% disp(l40)
% disp(l05)
% disp(l100) 
% disp(v03)
% disp(v0075)
% disp(v15)

%saccadeIndexNoFiltered = saccadeIndex;
save([outPath structname],'saccadeIndex', 'saccadeIndexAll','-v7.3')

%% use this to create a mat file and struct index containing a subset of the data
% with the given range of curvature, speed error, length, duration.
% For example, saccadeData_batch1 and saccadeData_batch1_Index are derived
% from the data in allSaccades3.mat and use the criteria below

% for i=1:length(tim)
%     dur(i,1)=tim{i}(end)-tim{i}(1);
% end
% sind=find(medcvn<.15 & err<.004 & saclength>4 & dur>=20);
% 
% eyeXY=sac(sind);
% spd=vel(sind);
% spdfit=velfit(sind);
% spdfiterr=err(sind);
% crv=medcvn(sind);
% len=saclength(sind);
% 
% save saccadeData_batch1 eyeXY spd spdfit spdfiterr crv len sind
% 
% saccadeIndex=struct;
% for i=1:length(sind)
%     saccadeIndex(i).saccade=sac{sind(i)};
%     saccadeIndex(i).videoNumber=vnum(sind(i));
%     saccadeIndex(i).trialNumber=trialnum(sind(i));
%     saccadeIndex(i).subjectNumber=subnum(sind(i));
%     saccadeIndex(i).elapsedtime=tim{sind(i)};
%     saccadeIndex(i).curvature=medcvn(sind(i));
%     saccadeIndex(i).length=saclength(sind(i));
%     saccadeIndex(i).direction=direction(sind(i));
%     saccadeIndex(i).speed=vel{sind(i)};
%     saccadeIndex(i).speedFit=velfit{sind(i)};
%     saccadeIndex(i).speedFitErr=err(sind(i));
% end
% 
% save saccadeData_batch1_Index saccadeIndex
