function [XY,sacind,eltimint,fixpos, XYnofilter]=offlineSaccadeDetect_noSmooth(eyedata,time)

%% thresholds
sthold=0.03;%0.15%.03; % speed threshold = 30 deg/sec
athold=0.02;%0.01;   % acceleration threshold
            % allSaccades3.mat and saccadeData_batch1 use .025 (25,000 deg/sec^2) 
            % Eyelink uses .008 (8,000 deg/sec^2)
       
%% get elapsed time
eltime=time-time(1);

% get interpolated elapsed time
eltimint=0:eltime(end);

eltime(1) = time(2)-time(1) - 1;


dx = diff(eltime);
indexes = find(~dx);

for i=1:length(indexes)
    eltime(indexes(i)) = eltime(indexes(i))-1;
end
   
eltimint = 1:length(time);

dx = diff(eltime);
indexes = find(~dx);

eltime(indexes) = []; %
eyedata(indexes,:) = [];

% interpolate eye position
warning('off','MATLAB:chckxy:IgnoreNaN');
x=spline(eltime,eyedata(:,1));
y=spline(eltime,eyedata(:,2));
xint=ppval(x,eltimint);
yint=ppval(y,eltimint);

 
% smooth the interpolated eye positions
% n=26;
% Wn=0.001;
% b=fir1(n,Wn);
% xfil=filter(b,1,[xint zeros(1,n/2)]);
% yfil=filter(b,1,[yint zeros(1,n/2)]);
% XY=[xfil(1+n/2:end)' yfil(1+n/2:end)'];
% XY = [xint' yint'];

X = sgolayfilt(xint,3,41);
Y = sgolayfilt(yint,3,41);
XY = [X' Y'];

% get velocity (time differential = 1 ms)
vel=[0 0;diff(XY)]; % interpolated data

% get speed
spd=sqrt(sum(vel.^2,2)); 

% remove blinks
spd(XY(:,1)<0)=0; % set missing data to 0
if any(spd==0)
   spd=removeBlinks2(spd,sthold);
end

% higher speed derivatives 
acc=[0 diff(spd')];

% get indices of above threshold values
spdind=spd>sthold;
accind=acc>athold;

wnsz=5; % 10 ms window
i=wnsz;
sacstart=0;
sacind=zeros(length(spdind),1);
while 1
   
    if i<length(spdind)-1
        i=i+1;
    else
        break
    end
    
    % The 2:6 change may help. Use if still getting small saccades; Raise
    % the acceleration threshold may help but not advisable
    if sacstart
        if  ~spdind(i-(wnsz-1)) %|| ...                   % if this datapoint is below speed threshold 
                           %(acc(i-(wnsz-1))<0 && acc(i-(wnsz-2))>0)    % or if acceleration passes above zero on the next sample;
                 
                % (acc(i-(wnsz-1))<0 && all(acc(i-(wnsz-(2:6))))>0)    % or if acceleration passes above zero on the next sample         
                 
            
            sacstart=0;                                 % signal end of saccade
            endind=i-(wnsz-1);
            sacind(startind:endind)=1;
            startind=[];
            continue
        else
            continue
        end
    elseif ~any(sacind(i-(wnsz):i))  &&...    % if the previous saccade didn't just end 
           all(spdind(i-(wnsz-1):i)) %&&...    % and if all datapoints in window are above speed threshold
         %  any(accind(i-(wnsz-1):i))          % and if the acceleration threshold is exceeded in the window
        sacstart=1;                           % signal start of saccade
        startind=i-(wnsz-1);
        endind=[];
        continue
    else
        continue
    end
end
%figure; plot(XY)
%hold on;plot(sacind*60,'r')
fixpos=XY(logical(1-sacind),:); % indices of fixations
shifts = (logical(1-sacind));
numberfixs = sum(abs(diff(logical(1-sacind))))/2;
fixinfo = [];
start = 1;
fixindex = 0;
startindex = 1;
X(X<0)=nan;
Y(Y<0)=nan;
X(X>100)=nan;
Y(Y>100)=nan;
for i=1:length(shifts)
    if shifts(i) == 1 && ~start       
       start = 1;
       startindex = i;
        
    else if ~shifts(i) && start
            start = 0;
            endindex = i-1;
            timing = endindex - startindex+1;
            if timing > 100 &&  timing<2500 
                fixindex = fixindex +1;                
                fixinfo(fixindex,1) = nanmean( X(startindex:endindex));
                fixinfo(fixindex,2) = nanmean( Y(startindex:endindex));
                fixinfo(fixindex,3) = timing;
            end
            
        end
    end
end
    
fixpos = fixinfo;

% In this version we do not smooth
XYnofilter = [xint' yint'];