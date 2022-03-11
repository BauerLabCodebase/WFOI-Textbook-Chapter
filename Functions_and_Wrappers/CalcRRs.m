function [R,Rs,seedradpix] = CalcRRs(refseeds,data,xform_isbrain)
nVy = size(data,1);
nVx = size(data,2);

mm=10;
mpp=mm/size(data,1);
seedradmm=0.25;
seedradpix=seedradmm/mpp;

if size(refseeds,1)>3
    for n=1:2:size(refseeds,1)-1
        if xform_isbrain(refseeds(n,2),refseeds(n,1))==1 && xform_isbrain(refseeds(n+1,2),refseeds(n+1,1))==1;  %%remove 129- on y coordinate for newer data sets
            SeedsUsed(n,:)=refseeds(n,:);
            SeedsUsed(n+1,:)=refseeds(n+1,:);
        else
            SeedsUsed(n,:)=[NaN, NaN];
            SeedsUsed(n+1,:)=[NaN, NaN];
        end
    end
else
    SeedsUsed= refseeds;
end
P=burnseeds(SeedsUsed,seedradpix,xform_isbrain);          % make a mask of seed locations
length_seeds=size(SeedsUsed,1);                           %reorder so bilateral connectivity values are on the off-diagonal
map=[(1:2:length_seeds-1) (2:2:length_seeds)];

data(isnan(data)) = 0;
data(isinf(data)) = 0;
data=reshape(data,nVy*nVx,[]);

disp('Calculating R Rs')
strace=P2strace(P,data,SeedsUsed);             % strace is each seeds trace resulting from averaging the pixels within a seed region
R= strace2R(strace,data , nVx, nVy);                 % R is the functional connectivity maps: normalize rows in time, dot product of those rows with strce

idx=find(isnan(SeedsUsed(:,1)));
R(:,:,idx)=0;


if size(refseeds,1)>3
    R_AP = makeRs(data ,strace);                     % R_AP is a matrix of the seed-to-seed fc values
    clear dataBandpass
    R_AP(idx, idx)=0;
    Rs = single(R_AP(map, map));
    clear strace R_AP
else
    Rs = [];
end

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             



