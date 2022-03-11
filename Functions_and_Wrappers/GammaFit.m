function [T,W,A,r,r2,hemoPred,objective_vals] = interSpeciesGammaFit(neural,hemo,t,mask)


%UNTITLED2 Summary of this function goes here
%   neural = 3D
%   hemo = 3D
% Output:
%   T = time to peak
%   W = FWHM
%   A = amplitude of gamma fit
%   r = correlation
%   r2 = goodness of fit
%   hemoPred = predicted hemodynamics


%Reshape
    OGsize=size(hemo);
    hemo=reshape(hemo,OGsize(1)*OGsize(2),[]);
    neural=reshape(neural,OGsize(1)*OGsize(2),[]);
    mask=reshape(mask,OGsize(1)*OGsize(2),[]);
    onlyBrain=find(mask==1);
    
    %Masking
    hemo=hemo(onlyBrain,:);
    neural=neural(onlyBrain,:);    
%Initialize
    hrfParam = nan(OGsize(1)*OGsize(2),3);
    hemoPred = nan(OGsize(1)*OGsize(2),OGsize(3));
    r = nan(1,OGsize(1)*OGsize(2));
    r2 = nan(1,OGsize(1)*OGsize(2));
    a=1;
    objective_vals=nan(1,OGsize(1)*OGsize(2));

% options = optimset('Display','iter','TolX',);

options=optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000,'TolX',1e-6,'TolF',1e-6);

for isbrain_ind=1:sum(mask)
            pixHemo = squeeze(hemo(isbrain_ind,:))';
            pixNeural = squeeze(neural(isbrain_ind,:))';
            he = HemodynamicsError(t,pixNeural,pixHemo);
            
%             worstErr = sum(pixHemo.^2);
%             options.TolFun = worstErr*0.01;
            
            fcn = @(param)he.fcn(param);
            
            
            [x,pixHrfParam,objective_val,exitflag,outputs] = evalc('fminsearchbnd(fcn,[1,1.5,0.25],[0.01,0.3,0.05],[2,3,1],options)'); %T, W, A -- guess, lower bound, upper bound

            pixelHrf = hrfGamma(t,pixHrfParam(1),pixHrfParam(2),pixHrfParam(3)); % %T, W, A
            pixHemoPred = conv(pixNeural,pixelHrf); 
            pixHemoPred = pixHemoPred(1:numel(pixNeural));


            hrfParam(onlyBrain(isbrain_ind),:) = pixHrfParam;
            hemoPred(onlyBrain(isbrain_ind),:) = pixHemoPred;
            r(onlyBrain(isbrain_ind)) = corr(pixHemoPred,pixHemo);%real(atanh(corr(pixHemoPred',pixHemo')));
            r2(onlyBrain(isbrain_ind))= 1-sumsqr(pixHemo-pixHemoPred)/sumsqr(pixHemo-mean(pixHemo));                           %1 - var(pixHemoPred - pixHemo)/var(pixHemo);
            objective_vals(onlyBrain(isbrain_ind))=objective_val;
end
    
%Reshape everything!
    r=reshape(r,OGsize(1),OGsize(2),[]);
    r2=reshape(r2,OGsize(1),OGsize(2),[]);
    objective_vals=reshape(objective_vals,OGsize(1),OGsize(2),[]);
    hemoPred=reshape(hemoPred,OGsize(1),OGsize(2),[]);
    hrfParam=reshape(hrfParam,OGsize(1),OGsize(2),[]);
 
%Extracting T,W,A
    T = hrfParam(:,:,1);
    W = hrfParam(:,:,2);
    A = hrfParam(:,:,3);

end