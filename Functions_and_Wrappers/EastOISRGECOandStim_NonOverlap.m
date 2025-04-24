function EastOISRGECOandStim_NonOverlap(ST,Stim)
%% AO controltomoboy
outV=5;     % output trigger height (volts)
Sr=10000;   % sampling rate

s = daq.createSession('ni');
s.Rate = Sr;

addAnalogOutputChannel(s, 'Dev1', 'ao0', 'Voltage');%Zyla #1,2
addAnalogOutputChannel(s, 'Dev1', 'ao2', 'Voltage');%M470
addAnalogOutputChannel(s, 'Dev1', 'ao3', 'Voltage');%M530
addAnalogOutputChannel(s, 'Dev1', 'ao4', 'Voltage');%TL625
addAnalogOutputChannel(s, 'Dev1', 'ao5', 'Voltage');%Stim
addAnalogOutputChannel(s, 'Dev1', 'ao6', 'Voltage');%pupil

numChannels = size(s.Channels,2);
stimChannel=5;%xw_220508 change numChannels to 5;


%% Control Matrix
DarkFrameTime=1;
LEDtime=round([90 3 90 6]); %blue, Green, Bright green,Red
readOutTime = 51;%minimum is 49
clearTime = 9;%29
SpF_Green = 100;%minimum is 48+49
SpF_Red = 100;%max(100,LEDtime(3)+2+readOutTime);
SpF_RGECO = LEDtime(3)+clearTime+readOutTime;
SpF_FAD = LEDtime(1)+clearTime+readOutTime;
SpF = [SpF_FAD SpF_Green SpF_RGECO SpF_Red ];

DataVec=zeros(sum(SpF),numChannels);
%laser on every other frame.Laser stim = 10hz, frame rate = 20hz;
CCDVec = zeros(sum(SpF),1);
CCDVec(1:SpF(1)-readOutTime,1) = outV;% FAD frame
CCDVec(SpF(1)+1:SpF(1)+SpF(2)-readOutTime,1) = outV; % green Frame
CCDVec(sum(SpF(1:2))+1:sum(SpF(1:2))+SpF(3)-readOutTime,1) = outV; % RGECO Frame
CCDVec(sum(SpF(1:3))+1:sum(SpF(1:3))+SpF(4)-readOutTime,1) = outV; % Red Frame


BlueVec = zeros(sum(SpF),1);
BlueVec(round((SpF(1)-readOutTime)/2-LEDtime(1)/2):round((SpF(1)-readOutTime)/2+LEDtime(1)/2),1) = outV;% blue OIS


GreenVec = zeros(sum(SpF),1);
GreenVec(SpF(1)+round((SpF(2)-readOutTime)/2-LEDtime(2)/2):SpF(1)+round((SpF(2)-readOutTime)/2+LEDtime(2)/2),1) = outV;% Green OIS
GreenVec(sum(SpF(1:2))+round((SpF(3)-readOutTime)/2-LEDtime(3)/2):sum(SpF(1:2))+round((SpF(3)-readOutTime)/2+LEDtime(3)/2),1) = outV;% Green bright

RedVec = zeros(sum(SpF),1);
RedVec(sum(SpF(1:3))+round((SpF(4)-readOutTime)/2-LEDtime(4)/2):sum(SpF(1:3))+round((SpF(4)-readOutTime)/2+LEDtime(4)/2),1) = outV;% RedOIS

DataVec(:,1) = CCDVec;
DataVec(:,2) = BlueVec;
DataVec(:,3) = GreenVec;
DataVec(:,4) = RedVec;

[time] = size(DataVec,1);

%% Stimulus Settings
switch Stim
    
    case 'On'           % Use for stimulating in a block design
        BL1=5;          % Baseline time (stimulus off)
        Dur=1;         % Stimulus duration time (stimulus on), doesnt matter as long as block length is longer than box' duration of being on
        BL2=30-(BL1+Dur);
        frames_block = round(BL1*Sr/(sum(SpF)))+round(Dur*Sr/(sum(SpF)))+round(BL2*Sr/(sum(SpF)));
        disp(['Frames for each block are ', num2str(frames_block)])%Xiaodan
        disp(['Frames for baseline of each block are ',num2str(round(BL1*Sr/(sum(SpF))))]);%Xiaodan
        
    case 'Off'              % No Stimuli
        BL1=60;
        Dur=0;
        BL2=0;
end

DataVecStim=DataVec;
DataVecStim(:,stimChannel)=outV;
DataVecNoStim=DataVec;
DataVecNoStim(:,stimChannel)=0.2;% xw 220508 from 0 to 0.2 since GalvoLaserOff is 0.2

%% Concatenate Everything Together
DataVecOff1=repmat(DataVecNoStim,round(BL1*Sr/sum(SpF)),1);
DataVecOn=repmat(DataVecStim,round(Dur*Sr/sum(SpF)),1);
DataVecOff2=repmat(DataVecNoStim,round(BL2*Sr/sum(SpF)),1);


DataVecTotal=cat(1, DataVecOff1, DataVecOn, DataVecOff2);
assignin('base','DataVecTotal', DataVecTotal);
Repeats=single(round(ST/(BL1+Dur+BL2)));
assignin('base','Repeats',Repeats);
DataVecTotal=repmat(DataVecTotal,Repeats,1);

DarkFrames=DataVec;
DarkFrames(:,2:numChannels)=0;
DarkFrames2=repmat(DarkFrames,round(DarkFrameTime*Sr/sum(SpF)),1);

DataVecTotal2 = cat(1,DarkFrames2,DataVecTotal); %ADAM
%DataVecTotal2(1,1)=outV; %for triggering first frame

numtrig=size(find(DataVecTotal2(:,1)==0),1)/readOutTime;

%% Begin Scanning

disp(['Set camera to record ',num2str(numtrig),' frames.'])
disp(['DarkFrame/Channel = ' num2str(DarkFrameTime*Sr./time)])
disp(['Frame rate of scan =  ',num2str(Sr./time),' fps.'])
disp('Press EN1TER when ready to start flashing sequence.')
disp('Playing encoding sequence.')
assignin('base','DataVecTotal2',DataVecTotal2);pause();pause(4)
queueOutputData(s, DataVecTotal2);
startForeground(s);


%% End scan

DataVec=zeros(Sr,size(s.Channels,2));
DataVec(:,1:size(s.Channels,2))=0;
queueOutputData(s, DataVec);
startForeground(s);
fclose('all');
beep on; beep
end