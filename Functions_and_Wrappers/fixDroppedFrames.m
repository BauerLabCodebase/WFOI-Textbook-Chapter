function [dfIndRaw, dfIndFixed, fixedRaw, fixedRawTime] = fixDroppedFrames(runInfo,raw, rawTime, fixedData)

% load rawQCData
currMouseName = [runInfo.recDate '-' runInfo.mouseName '-' runInfo.session...
    '' num2str(runInfo.run)];
if fixedData
    rawQCDatLoc = [runInfo.saveFolder filesep currMouseName '-rawQCData_dropCorr.mat'];
else
    rawQCDatLoc = [runInfo.saveFolder filesep currMouseName '-rawQCData.mat'];
end
rawQCDat = load(rawQCDatLoc);

% check if dropped frames exist
instChange = rawQCDat.diffF2F;
chData = rawQCDat.rawSpatialAvg;

[instChangeVal,instChangeInd]=maxk(instChange,25); % look at 25 largest maximums in inst change

dfThreshold = 1e7; % threshold set for minimum inst change required to detect dropped frame
dfInd = []; % store indicies of dropped frames


for i=1:length(instChangeVal)
    if instChangeVal(i) > dfThreshold
        dfInd = [dfInd instChangeInd(i)];
    end
end

dfInd = sort(dfInd);
dfIndRaw = dfInd - 1; %subtract one to account for shift due to inst change calculations

if isempty(dfIndRaw) % if no df indicies are found, end function since no correction is needed
    dfIndFixed = dfIndRaw;
    fixedRaw = [];
    fixedRawTime = []; 
    return;
else
    % crop dropped frames within light level data to get df indicies in fixed data
    dfIndFixed = dfIndRaw;
    for i=1:length(dfIndFixed)
        dropInd = dfIndFixed(i);
        chData(dropInd,:) = []; % removed dropped frames
        dfIndFixed(i+1:end) = dfIndFixed(i+1:end) - 1;
    end
    dfIndFixed = unique(dfIndFixed); % only keep unique indicies
    
    % use light level data to get permutation matrix
    avgB = zeros(length(dfIndFixed),4); % store avg before dropped frame
    avgA = zeros(length(dfIndFixed),4); % store avg after dropped frame
    numFrames = 11; % how many frames to look at before and after dropped frame
    permMat = zeros(length(dfIndFixed),4); % initialize permutation matrix
%     permMat = [2 3 4 1];
    
    for j=1:length(dfIndFixed)
        % check to see if there are enough data points before and 
        % after the dropped frame, and, if not, crop the excess data
        if (dfIndFixed(j)+numFrames) > length(chData(:,1))
            raw(:,:,:,dfIndFixed(j):end) = [];
            dfIndFixed(j) = []; %possible paradox
        elseif (dfIndFixed(j)-numFrames) <= 1
            disp('Error: early frame drop detected!');
            pause;
        else
            for i=1:4 % four channels
                avgB(j,i) = mean(chData(dfIndFixed(j)-numFrames:dfIndFixed(j)-1,i));
                avgA(j,i) = mean(chData(dfIndFixed(j)+1:dfIndFixed(j)+numFrames,i));
            end

            % determine which channel does original ch1 switch to
            diffMat = abs(abs(avgA(j,:)-avgB(j,1)));
            [minDiff,newCh] = min(abs(avgA(j,:)-avgB(j,1)));
            
            % check to see if channles have similar differences
            % this is a check to see if two traces are too close to
            % distinguish properly
            minDelta = diffMat - minDiff;
            simVals = minDelta(minDelta>0 & minDelta<0.025e3); % similar values
            if ~isempty(simVals)
                simCh = [];
                for jj = 1:length(simVals)
                    simCh = [simCh find(minDelta==simVals(jj))];
                end
                disp(['For frame drop #' num2str(j) ', channel(s) ' num2str(simCh) ' counts are too close to Ch '...
                    num2str(newCh) ' for proper distinciton.']);
                qcimg = [runInfo.saveRawQCFig '.png'];
                qcimg = imread(qcimg);
                figure(1);
                set(figure(1),'Position',[100 100 1400 900]);
                zoom on;
                openQCimg = image(qcimg);
                userPrompt = 'Please look at Figure 1 and enter the channel that joins to Ch1 after frame drop: ';
                userInput = input(userPrompt);
                close(figure(1));
                perm = [userInput userInput+1 userInput+2 userInput+3];
                permInd = find(perm>4);
                perm(permInd) = perm(permInd)-4;
                permMat(j,:) = perm; % store permutation values in matrix               
            else
                perm = [newCh newCh+1 newCh+2 newCh+3];
                permInd = find(perm>4);
                perm(permInd) = perm(permInd)-4;
                permMat(j,:) = perm; % store permutation values in matrix 
            end
        end
    end 
end
% perform frame drop correcting on actual image data
% crop dropped frames from raw data
dfIndFixed = dfIndRaw;
for i=1:length(dfIndFixed)
    dropInd = dfIndFixed(i);
    if dropInd < length(raw) %if it is not the last index
         %for all channels
            raw(:,:,:,dropInd) = [];      
    end
    dfIndFixed(i+1:end) = dfIndFixed(i+1:end) - 1;
end

dfIndFixed = unique(dfIndFixed);

% perform correction for channel switch
fixedRaw = raw;
for j=1:length(dfIndFixed) 
    if dfIndFixed(j) < length(raw) %if it is nt the last index
         for i=1:4
             fixedRaw(:,:,i,dfIndFixed(j):end) = raw(:,:,permMat(1,i),dfIndFixed(j):end);
         end
    end
    raw = fixedRaw;
end

% adjust time scale to account for cropped frames
fixedRawTime = rawTime;
fixedRawTime(:,(length(fixedRaw(:,:,1,:))+1):end)= []; 
end