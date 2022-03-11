files = dir('Y:\RAD\**\*Landmarks*.mat');
%this script allows you to obtain a white light image of the brain, with
%boundaries marked out by grayscale and landmarks indicated with dots. It
%saves all overlays as a .png file with the same location and name as the
%original "LandmarksAndMask.mat" file. 


for i=1:length(files)
    fileboy=[files(i).folder '\' files(i).name];
    load(fileboy)
    savename=fileboy(1:end-3);
    
    imagesc(~isbrain.*rgb2gray(WL))
    colormap('Gray')
    hold on;
    
    imagesc(WL,'alphadata',isbrain);hold on;
    plot(I.tent(1)*128,I.tent(2)*128,'o')
    plot(I.bregma(1)*128,I.bregma(2)*128,'o')
    plot(I.OF(1)*128,I.OF(2)*128,'o')
    saveas(gcf,[savename 'png'])
    close all;
end

