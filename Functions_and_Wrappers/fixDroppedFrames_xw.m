function [dfIndRaw, dfIndFixed, fixedRaw] = fixDroppedFrames_xw(runInfo,raw,fixedData,subfile)

% load rawQCData
currMouseName = [runInfo.recDate '-' runInfo.mouseName '-' runInfo.session...
    '' num2str(runInfo.run)];
if fixedData
    rawQCDatLoc = [runInfo.saveFolder filesep currMouseName subfile '-rawQCData_dropCorr.mat'];
else
    rawQCDatLoc = [runInfo.saveFolder filesep currMouseName subfile '-rawQCData.mat'];
end
rawQCDat = load(rawQCDatLoc);

% check if dropped frames exist
instChange = rawQCDat.diffF2F;
chData = rawQCDat.rawSpatialAvg;

[instChangeVal,instChangeInd]=maxk(instChange,25); % look at 25 largest maximums in inst change

dfThreshold = 1e9; % xw_threshold set for minimum inst change required to detect dropped frame
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
    return;
else
    % crop dropped frames within light level data to get df indicies in fixed data
    dfIndFixed = dfIndRaw;
    dfIndFixed = unique(dfIndFixed); % only keep unique indicies
    
    % use light level data to get permutation matrix
    avgB = zeros(length(dfIndFixed),3); % store avg before dropped frame
    numFrames = 11; % how many frames to look at before and after dropped frame

    
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
            for i=1:3 % three channels
                avgB(j,i) = mean(chData(dfIndFixed(j)-numFrames:dfIndFixed(j)-1,i));
                avgD(j,i) = chData(dfIndFixed(j),i);
            end

            % determine which channel does original ch1 switch to
            [~,newCh_1] = min(abs(avgD(j,1)-avgB(j,:)));
            % determine which channel does original ch2 switch to
            [~,newCh_2] = min(abs(avgD(j,2)-avgB(j,:)));
            % determine which channel does original ch3 switch to
            [~,newCh_3] = min(abs(avgD(j,3)-avgB(j,:)));
            
        end
    end 
end
% perform frame drop correcting on actual image data
% make droped frames same as the frames before it.
% perform correction for channel switch
fixedRaw = raw;
dfIndFixed = dfIndRaw;
dfIndFixed = unique(dfIndFixed);
fixedRaw = reshape(fixedRaw,128,128,[]);
raw = reshape(fixedRaw,128,128,[]);
for i=1:length(dfIndFixed)
    dropInd = dfIndFixed(i);
    if dropInd < length(raw/runInfo.numCh) %if it is not the last index
         %for all channels
         if newCh_1 == 1 && newCh_2 == 3 && newCh_3 == 1
         fixedRaw(:,:,(dropInd-1)*runInfo.numCh+2) = raw(:,:,(dropInd-1)*runInfo.numCh-1);
         fixedRaw(:,:,(dropInd-1)*runInfo.numCh+3:end) = raw(:,:,(dropInd-1)*runInfo.numCh+2:end-1);
         elseif newCh_1 == 2 && newCh_2 == 3 && newCh_3 == 1
         fixedRaw(:,:,(dropInd-1)*runInfo.numCh+1) = raw(:,:,(dropInd-1)*runInfo.numCh-2);
         fixedRaw(:,:,(dropInd-1)*runInfo.numCh+2:end) = raw(:,:,(dropInd-1)*runInfo.numCh+1:end-1);             
         elseif newCh_1 == 1 && newCh_2 == 2 && newCh_3 == 1
         fixedRaw(:,:,(dropInd-1)*runInfo.numCh+3) = raw(:,:,(dropInd-1)*runInfo.numCh-3);
         fixedRaw(:,:,(dropInd-1)*runInfo.numCh+4:end) = raw(:,:,(dropInd-1)*runInfo.numCh+3:end-1);              
         end
    end
    raw = fixedRaw;
end
fixedRaw = reshape(fixedRaw,128,128,runInfo.numCh,[]);
end