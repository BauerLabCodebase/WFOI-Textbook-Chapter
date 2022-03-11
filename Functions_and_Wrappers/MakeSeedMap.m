function seedMap=MakeSeedMap(numPixX,numPixY,seedCenter,seedRadius)
    seedNum = size(seedCenter,1);
    seedMap = false(numPixY,numPixX,seedNum);
    for seedInd = 1:seedNum
        seedCoor = circleCoor(seedCenter(seedInd,:),seedRadius);
        seedCoor = matCoor2Ind(seedCoor,[numPixY numPixX]);
        seedBool = false(numPixY,numPixX); seedBool(seedCoor) = true;
        seedMap(:,:,seedInd) = seedBool;
    end

end