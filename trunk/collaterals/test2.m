classdef test2
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        c
    end
    
    methods
        function obj = test2(input)
            obj.c = obj.calc(input)
        end
        
    end
    methods(Static)
        function  a= calc(d)
            a =d^2;
        end
    end
    
end

