clear all
close all
k=0;
sac={};
tim={};
vnum=[];
trialnum=[];
subnum=[];
limit1 = 25000;
limit2 = 5010;

%type of DB %1 BioEye 2 Gaze
DBtype = 1;

files_to_read = 1;
switch DBtype
    
    case 0
        etPath='~/Desktop/AckermannArchive/Video clip eyetracking data';
        load 'video number lookup.mat'
        load baselineSubjectIndex
        screenDist=45;        % distance in cm
        screenSize=[60 34];   % [horz vert] screen dimensions in cm
        screenRes=[2560 1440];% screen dims in pixs
        files_to_read = max(videoNumbers);
        
    case 1
        % BioEye files
        direct = '~/Desktop/GazeContigentDisplay/BioEye2015_DevSets/RAN_1year_dv/';
        direct2 ='~/Desktop/GazeContigentDisplay/BioEye2015_DevSets/RAN_30min_dv/';
        files = dir([direct '*.txt']);
        files2 = dir([direct2 '*.txt']);
        files3 = [files ; files2];
        %filename = '~/Desktop/GazeContigentDisplay/BioEye2015_DevSets/RAN_1year_dv/ID_006_1.txt';
        files_to_read = length(files3);
        samplerate = 250;
        screenDist=55;        % distance in cm474×297 mm, and screen resolution of 1680×1050 pixels.
        screenSize=[47.4 29.7];   % [horz vert] screen dimensions in cm
        screenRes=[1680 1050];% screen dims in pixs
    case 2
        % Gaze files Michael
        direct = '~/Desktop/GazeContigentDisplay/gaze/natural_movies_gaze/';
        filename = [direct 'ALK_ducks_children.coord'];
        files3 = dir(filename);
        samplerate = 250;         
        screenDist=45;        % distance in cm
        screenSize=[40 22.5];   % [horz vert] screen dimensions in cm
        screenRes=[1280 720];% screen dims in pixs
        files_to_read = length(files3);
        
end

saccadeIndex=struct;
for fil_index =1:files_to_read
    
    %disp(files3(fil_index).name);
    if DBtype > 0
        % open the file
        if fil_index > 74 && DBtype==1
            direct = direct2;
        end
        filename = [direct files3(fil_index).name];
        fid = fopen(filename);
        finfo = dir(filename);
        fsize = finfo.bytes;
        D = textscan(fid, '%s','delimiter', '\n');
        
%         screenDist=45;        % distance in cm
%         screenSize=[40 22.5];   % [horz vert] screen dimensions in cm
%         screenRes=[1280 720];% screen dims in pixs
        
        switch DBtype
            case 1
                for i=2:length(D{1})
                    values(i-1,:) = sscanf(D{1}{i}, '%g %g %g %g %g %g',6);
                end
                 subjectNum = str2num(files3(fil_index).name(4:6));
                 numTrial = str2num(files3(fil_index).name(8));
            case 2
                for i=3:length(D{1})
                    values(i-2,:) = sscanf(D{1}{i}, '%g %g %g %g',4);
                    values(i-2,2:3) = pix2deg([values(i-2,2)' values(i-2,3)'],screenDist,screenSize,screenRes);
                    
                end
                 subjectNum = files3(fil_index).name(1:3);
        end
        % read each set of measurements
        
        % close the file
        fclose(fid);
        
        ed=[values(:,2) values(:,3)];
        time = values(:,1)';
        
        
        switch DBtype
            case 1
                [edint{1},saccind{1},timint{1},fixpos, edintNoFilter{1}]=offlineSaccadeDetectModified250(ed,time);
            case 2
                [edint{1},saccind{1},timint{1},fixpos]=offlineSaccadeDetectModified400(ed,time);
        end
       
        
    else
        
        trials4Vid=find(videoNumbers==fil_index);
        
        % get the subject number for each trial
        subs4Vid=subjectIndex(trials4Vid); % see subjectList for sub IDs
        
        % get all eye data for this clip
        for i = 1:length(trials4Vid)
            load([etPath filesep eyetrackFiles{trials4Vid(i)}]);
            temp.x = eyetrackRecord.x;
            temp.y = eyetrackRecord.y;
            temp.t = eyetrackRecord.t;
            temp.missing = eyetrackRecord.missing;
            etData(i) = temp;
        end
        
        % loop through each trial
        edint=cell(length(etData),1);
        saccind=cell(length(etData),1);
        timint=cell(length(etData),1);
        trialind=zeros(length(etData),1);
        subind=zeros(length(etData),1);
        for i=1:length(etData)
            if length(etData(i).t)>19000 % make sure there's enough data
                ed=[etData(i).x' etData(i).y'];
                time=etData(i).t;
                % convert to degrees
                eddeg=pix2deg(ed,screenDist,screenSize,screenRes);
                % index saccades
                [edint{i},saccind{i},timint{i},fixpos]=offlineSaccadeDetect4(eddeg,time); % original
                % index trial
                trialind(i,1)=trials4Vid(i);
                % index subject
                subind(i,1)=subs4Vid(i);
            end
        end
        
    end
    
   
   
    % plots
    % original data
    %     figure;set(gcf,'color',[1 1 1])
    %     plot(ed(:,1),'r');hold on
    %     plot(ed(:,2),'g')
    % original data w/o blinks
    % interpolated data and saccade index
    %     figure;set(gcf,'color',[1 1 1])
    %     plot(edint{1}(:,1),'r');hold on
    %     plot(edint{1}(:,2),'g')
    %     V=axis;
    %     plot(V(end)*saccind{1},'k')
    %     xlabel('time (ms)')
    %     ylabel('position (pix)')
    %     legend({'horz' 'vert' 'saccade index'})
    
    % figure;plot(values(:,2),'b');hold on
    % plot(saccind{1},'k');
    
    % if DBtype ==1
    % for i=250:250:250*round(limit/250)
    %     line([i,i],[-20,15],'Color',[1 0 0]);
    % end
    % end
    % extract the saccades
  
    % extract the saccades
    for sub=1:length(edint)
        insac=false;
        for i=1:length(saccind{sub})
            if saccind{sub}(i)==1 && ~insac
                k=k+1;
                sac{k}=edint{sub}(i,:); % saccade coordinates
                sacNoFilter{k}=edintNoFilter{sub}(max(1,i-4):i,:); % saccade coordinates
                tim{k}=timint{sub}(i);  % get clock time of the saccade
                vnum(k)=fil_index;         % the clip number                               
               
                if DBtype<2
                    subnum(k)=subjectNum; %subind(sub);%subjectNum;
                    trialnum(k) = numTrial;
                else
                    subnum{k} = subjectNum;
                    trialnum(k)=trialind(sub);  % the trial it's derived from (1:2681)
                end
                insac=true;
                continue
            elseif saccind{sub}(i)==1 && insac
                sac{k}(end+1,:)=edint{sub}(i,:);
                sacNoFilter{k}(end+1,:)=edintNoFilter{sub}(i,:);
                tim{k}(end+1,1)=timint{sub}(i);
                continue
            elseif saccind{sub}(i)==0 && insac
                insac=false;
                continue
            elseif saccind{sub}(i)==0 && ~insac
                continue
            end
        end
    end 
    clear etData
   
    disp(length(sac));
end
   

    disp('curvature');
    %curvature
    medcvn=zeros(length(sac),1);
    for i=1:length(sac)
        X=sac{i}(:,1);
        Y=sac{i}(:,2);
%         if length(X)>9
%             X2 = interp(X,1000/samplerate);
%             Y2 = interp(Y, 1000/samplerate);
%         else
%             X2 = resample(X,1000/samplerate,1);
%             Y2 = resample(Y,1000/samplerate,1);
%         end
        [cv cvn]=curvature(X,Y);
        medcvn(i)=median(cvn);
    end
    
    
    disp('length');
    % get length
    saclength=zeros(length(sac),1);
    for i=1:length(sac)
        X=sac{i}(:,1);
        Y=sac{i}(:,2);
        saclength(i)=norm([X(end) Y(end)]-[X(1) Y(1)]);
    end        
    
    disp('direction');
    % get direction
    direction=zeros(length(sac),1);
    for sacind=1:length(sac)
        XY=sac{sacind};
        endpoint=XY(end,:)-XY(1,:);
        ang=acos(endpoint(1)/norm(endpoint));
        if endpoint(2)<0
            ang=2*pi-ang;
        end
        direction(sacind)=ang*180/pi;
    end
        
    disp('velocity');
    % get velocity of each saccade and fit exponential velocity fn
    vel=cell(length(sac),1);
    velfit=cell(length(sac),1);
    params=cell(length(sac),1);
    err=zeros(length(sac),1);
   
    for sacind=1:length(sac)
        XY=sac{sacind};
        if length(XY(:,1))>9
            X2 = interp(XY(:,1),1000/samplerate);
            Y2 = interp(XY(:,2), 1000/samplerate);
        else
            X2 = resample(XY(:,1),1000/samplerate,1);
            Y2 = resample(XY(:,2),1000/samplerate,1);
        end
        XY = [X2 Y2];
        vel{sacind}=sqrt(sum(diff(XY).^2,2));
        [velfit{sacind} params{sacind} err(sacind)]=fitVelocity(vel{sacind},4);
    end
    
    disp('Indexes');
    saccadeIndex = {};
    j = 1;
    for i=1:length(sac)          %
        %         if DBtype == 1
        %
        %             %saccadeIndex(i).trialNumber=trialnum(i);
        %         end
%         if i==95
%            disp(1) 
%         end
        if length(sac{i})>15 && saclength(i) <60 && saclength(i)>1 && length(sac{i})<150 && vel{i}(end)<0.3 &&  vel{i}(round(length(sac{i})/4)) > 0.15*max(vel{i})
            if DBtype ==1
                saccadeIndex(j).elapsedtime=tim{i};
                saccadeIndex(j).trialNumber=vnum(i);
            else
                saccadeIndex(j).elapsedtime=unique(round(tim{i}/1000));
            end
            if DBtype == 0
                saccadeIndex(j).videoNumber=vnum(i);
                saccadeIndex(j).trialNumber=trialnum(i);
            end
            saccadeIndex(j).saccade=sac{i};
            
            saccadeIndex(j).subNumber=subnum(i);
            saccadeIndex(j).saccadeNoFilter=sacNoFilter{i};
            saccadeIndex(j).curvature=medcvn(i);
            saccadeIndex(j).length=saclength(i);
            saccadeIndex(j).trialNumber=trialnum(i);
            saccadeIndex(j).direction=direction(i);
            saccadeIndex(j).velocity.velocity=vel{i};
            saccadeIndex(j).velocity.velocityFit=velfit{i};
            saccadeIndex(j).velocity.params=params{i};
            saccadeIndex(j).velocity.err=err(i);
            
            if DBtype==1
                saccadeIndex(j).duration= length(saccadeIndex(j).elapsedtime)*(1000/samplerate); 
            else
                saccadeIndex(j).duration= length(saccadeIndex(j).elapsedtime)*(1000/samplerate); 
            end
            j = j+1;
        end
    end
    
 
    
 save('SaccadeIndexBioEye.mat', 'saccadeIndex');
%     nu = 1;
%   figure
%  plot (values(:,2))
%  hold on
%  for i=1:length(saccadeIndex)
%  plot(saccadeIndex(i).elapsedtime,values(saccadeIndex(i).elapsedtime,2), 'r','LineWidth',7);
%  hold on
%   plot(saccadeIndex(i).elapsedtime,values(saccadeIndex(i).elapsedtime,1), 'b','LineWidth',7);
% 
%  end


% % Old plot type
% plot(ed(:,1));
% hold on;
% for i=1:length(saccadeIndex)
% plot(saccadeIndex(i).elapsedtime,ed((saccadeIndex(i).elapsedtime),1), 'r','LineWidth',6);
% end
% xlabel('Time (s)','FontSize', 18);
% ylabel('X trace','FontSize', 18);
% 
% set(gca,'FontSize', 18, 'XTick', 0:250:ed(end,1),  'XTickLabel', 0:1:ed(end,1));


%    save(['saccades' num2str(DBtype) '.mat'], 'saccadeIndex', 'sac', 'medcvn', 'err' , 'saclength', 'direction', 'velfit', 'vel');
%     cvrange=[.21 .30]%[.17 .3];
%     cvind=find(medcvn>=cvrange(1) & medcvn<cvrange(2))
%     figure;set(gcf,'color',[1 1 1])
%     for i=1:length(cvind)
%         x=sac{cvind(i)}(:,2)-sac{cvind(i)}(1,2);
%         y=sac{cvind(i)}(:,1)-sac{cvind(i)}(1,1);
%         plot(x,y,'o-')
%         set(gca,'xlim',[-40 40],'ylim',[-40 40],...
%             'xtick',[-40:5:40],'ytick',[-40:5:40])
%         xlabel('horz (deg)')
%         ylabel('vert (deg)')
%         text(-35,35,['curvature =' num2str(medcvn(cvind(i)))])
%         %     set(gca,'xlim',[-500 500],'ylim',[-500 500],...
%         %             'xtick',[-500:100:500],'ytick',[-500:100:500])
%         %     xlabel('horz (pix)')
%         %     ylabel('vert (pix)')
%         text(-450,450,['curvature =' num2str(medcvn(cvind(i)))])
%         axis square
%         grid on
%         pause
%     end
%     
%     %% plot saccade and velocity fit
%     %load allSaccades3
%     %load allSaccades3_velocity
%     sind=find(medcvn<.45 & err<.004); % It was 0.045 before
%     %sind=[5221 7593 7594 10073 11302 11670 15701 15752 16443]; % saccades that change direction
%     figure;set(gcf,'color',[1 1 1])                            % indices from allSaccades2
%     for i=1:length(sind)
%         %plot saccade
%         subplot(1,2,1)
%         x=sac{sind(i)}(:,2)-sac{sind(i)}(1,2);
%         y=sac{sind(i)}(:,1)-sac{sind(i)}(1,1);
%         plot(x,y,'o-')
%         %     set(gca,'xlim',[-700 700],'ylim',[-500 500],...
%         %             'xtick',[-700:100:700],'ytick',[-500:100:500])
%         set(gca,'xlim',[-20 20],'ylim',[-20 20],...
%             'xtick',[-20:2:20],'ytick',[-20:2:20])
%         axis square
%         grid on
%         xlabel('horz (deg)')
%         ylabel('vert (deg)')
%         text(-18,18,['curvature =' num2str(medcvn(sind(i)))])
%         
%         % plot velocity and fit
%         subplot(1,2,2)
%         plot(vel{sind(i)},'bo-','markersize',12);hold on
%         plot(velfit{sind(i)},'r','linewidth',2);hold off
%         set(gca,'ylim',[0 1],'ytick',0:.1:1)
%         xlabel('t (ms)')
%         ylabel('speed (deg/ms)')
%         title(['err = ' num2str(err(sind(i))) '   saccade # ' num2str(sind(i))])
%         pause
%     end
%     



