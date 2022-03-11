classdef HemodynamicsError
    %HemodynamicsError calculates the error difference between predicted
    %and actual hemodynamics
    %   Detailed explanation goes here
    
    properties
        NeuralActivity
        ActualHemodynamics
        Time
    end
    
    methods
        function obj = HemodynamicsError(time,neuralActivity,actualHemodynamics)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            
            if size(neuralActivity,1) ~= 1
                neuralActivity = neuralActivity';
            end
            if size(actualHemodynamics,1) ~= 1
                actualHemodynamics = actualHemodynamics';
            end
            if size(time,1) ~= 1
                time = time';
            end
            
            obj.NeuralActivity = neuralActivity;
            obj.ActualHemodynamics = actualHemodynamics;
            obj.Time = time;
        end
        
        function diffVal = fcn(obj,param)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            T = param(1); W = param(2); A = param(3);
            
            convData = conv(obj.NeuralActivity, hrfGamma(obj.Time,T,W,A));
            convData = convData(1:numel(obj.NeuralActivity));
            
            diffVal = sum((convData - obj.ActualHemodynamics).^2);
        end
    end
end

