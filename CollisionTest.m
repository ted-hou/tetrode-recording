classdef CollisionTest < handle
    properties
        Trials = [],
        Timestamps = [],
        StimOn = [],
        StimOff = [],
        TrainOn = [],
        TrainOff = [],
        Window = [0, 0],
        ElectrodesInfo = struct([])
        Filename = struct('NSx', '', 'TR', '', 'PTR', '')
    end

    methods
        function obj = CollisionTest()
            files = uipickFiles('Prompt', 'Select raw data folders.', 'Type')
        end
    end

    methods (Static)
        function obj = staticMethodName(args)
            
        end
    end
end