%% load in Atlas and Network Assignments
load('AtlasandIsbrain.mat', 'AtlasSeedsAllNetworksLR', 'NetworksLR', 'CmapAllNetworks')

%% load in movie and brain mask
load('A:\JPC\Old Lis\Old_data_run\201124\201124-G52M5baseline-fc3-dataFluor.mat')


%% Calculate whole cortex FC matrix
xform_dataflour_GSR=gsr(xform_datafluorCorr,xform_isbrain);
xform_dataflour_GSR=reshape(xform_dataflour_GSR,128*128,[]);
TempR=corr(xform_dataflour_GSR',xform_dataflour_GSR'); 

%% Symmetrize Brain Mask 
symisbrainall=zeros(128);
symisbrainall(:,1:64)=xform_isbrain(:,1:64).*fliplr(xform_isbrain(:,65:128)); %this multiplies both sides of the brain to get the intersection of both hemispheres. This is then reflected
symisbrainall(:,65:128)=fliplr(symisbrainall(:,1:64));
symisbrainall=uint8(symisbrainall);
symisbrainall=single(symisbrainall);
Brainidx=find(symisbrainall==1); %find the indices that are in this new mask

%% Crop Whole cortex FC matrix
TempR=TempR(Brainidx,Brainidx); %This makes the matrix much smaller!

%% Remove networks containing only a few pixels (we can discuss)
Atlas=AtlasSeedsAllNetworksLR;
BadNets=[23,24,4,1,25,21,13,12,8,3,7,22];
BadNets=[BadNets; BadNets+24];

for p=1:numel(BadNets)
    idx=Atlas==BadNets(p);
    Atlas(idx)=0;
end

%% Reorder indices of brain mask to go along rows instead of columns (only affects final look of whole cortex FC matrix)
APPixels=[];
a=0;
for y=1:128
    for x=1:128
        if symisbrainall(y,x) %only if it is in the brain mask
            a=a+1;
            APPixels(a)=sub2ind([128, 128], y, x);
        end
    end
end
APPixels=APPixels'; %Anterior-Posterior?
%% Make a mapping vector between the indicies in the brain mask and the rows of the whole cortex FC matrix

%Mapping between verticle unwrapping (Brainidx) and horizontal unwrapping
%(APPixelsidx)
APPixelsidx=[];
for n=1:size(APPixels,1)
    APPixelsidx(n)=find(Brainidx==APPixels(n),1, 'first'); 
end
APPixelsidx=APPixelsidx';

%% Reorder whole cortex correlation matrix
TempR=TempR(APPixelsidx,APPixelsidx);

%% find rows of correlation matrix that correspond to network parcels and ignore unassigned brain pixels
[Nets,I] = sort(Atlas(APPixels),'ascend'); % sort by network in Left then Right
I(Nets == 0) = [];
GoodNets=unique(Nets(Nets~=0))';
Nets(Nets == 0) = [];

%% How can we order in another way? For instance, 1,11,2,12..etc

target_order=[1,25,2,26,3,27,4,28,5,29,6,30,7,31,8,32,9,33,10,34,11,35,12,36,13,37,14,38,15,39,16,40,17,41,18,42,19,43,20,44,21,45,22,46,23,47,24,48];

I=[];
Nets=[];
NetworksLR_new=NetworksLR;

for i=1:numel(target_order)
   
        I=[I; find(Atlas(APPixels)==target_order(i))];
        Nets(find(Atlas(APPixels)==target_order(i)))= target_order(i)  ;
        NetworksLR_new(i,:)= NetworksLR(  find([NetworksLR{:,2}]==target_order(i)),  :  );
end
Nets=Nets';
NetworksLR=NetworksLR_new;
GoodNets=unique(Nets(Nets~=0))'; %this automatically sorts the indeices in ascending order, which is what we want for the earlier case. We want to preseve the order seen in target_order though...

%
BadNets=setdiff(target_order,GoodNets);

tmp=[];
for i=1:numel(BadNets)
    tmp=[tmp find(target_order==BadNets(i))];
end
GoodNets=target_order;
GoodNets(tmp)=[];


Nets(Nets == 0) = [];

%% Reorder according to network assignment
ReorderedFCMat=TempR(I,I);

%% Plot

clear key Cmap
%initializing variables we will need
key(:,1)=1:size(ReorderedFCMat,1);
key(:,2)=Nets;
CmapAllNetworksLR=[CmapAllNetworks; CmapAllNetworks];


%organizing based on the ORDER 
order=cell2mat(NetworksLR(:,2)); %this is all of the network assignments and the order you want them in?

OrderedNetworks(order,1)=join([NetworksLR(:,1), NetworksLR(:,3)]); %combining the network with the 'Left' or 'Right'. 
OrderedCmapAllNetworks(order,:)=CmapAllNetworksLR; %change the order of the colors 

GoodNetsV=[GoodNets; GoodNets];%?
GoodNetsH=[GoodNets; GoodNets]; %?
OrderedNetworksH=OrderedNetworks; %?

Cmap=OrderedCmapAllNetworks;

buffer=200;
lims=[-0.7 0.7];

figure; 
Matrix_Org3(ReorderedFCMat, key, GoodNetsV, GoodNetsH, OrderedNetworks, OrderedNetworks, buffer, lims, Cmap,'jet')
