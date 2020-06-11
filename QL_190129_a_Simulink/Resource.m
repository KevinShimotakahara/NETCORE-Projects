% This is a generic Resource object that
% also keeps track of statistics on the number of busy resources

classdef Resource < handle
    properties
        Busy %Integer
        NumberOfUnits %Integer
        NumBusy %CTStat
    end
    
    methods
        
        function obj = Resource(varargin)
            % Executes when resource object is created to initialize variables
            % and add number of busy Resources statistic to TheCTStats collection
            global MyGV
            switch nargin
                case 0
                    obj.Busy = 0;
                    obj.NumberOfUnits = 0;
                    obj.NumBusy = CTStat;
                    MyGV.TheCTStats = [MyGV.TheCTStats obj.NumBusy];
                otherwise
                    error('Cannot set property values')
            end
        end
        function result = Seize(obj, Units)
            % Seize Units of resource then updates statistics
            % Returns False and does not seize if not enough resources available;
            % otherwise returns True
            
            diff = obj.NumberOfUnits - Units - obj.Busy;
            % If diff is nonnegative, then there are enough resources to seize
            if diff >= 0
                obj.Busy = obj.Busy + Units;
                obj.NumBusy.Record(double(obj.Busy));
                result = true;
            else
                result = false;
            end
        end
        
        function result = Free(obj, Units)
            % Free Units of resource then updates statistics
            % Returns False and does not free if attempting to free more
            % resources than available; otherwise returns True
            % CORRECTED 2/1/2010
            
            diff = obj.Busy - Units;
            % If diff is negative, then trying to free too many resources
            if diff < 0
                result = false;
            else
                obj.Busy = obj.Busy - Units;
                obj.NumBusy.Record(double(obj.Busy));
                result = true;
            end
        end
        function result = Mean(obj)
            % Return time-average number of busy resources up to current time
            result = obj.NumBusy.Mean;
        end
        
        function SetUnits(obj, Units)
            % Set the capacity of the resource (number of identical units)
            
            obj.NumberOfUnits = Units;
        end
        
    end
end



