classdef constantDMdl < diffusion_function
    
    properties
        nmax
        rres
        rlim
        sigma
        finalConc
        D
    end
    
    methods
        function obj = constantDMdl(options)
            obj@diffusion_function(options);
        end
        function obj = calcCurve(obj,p)
            obj = obj.updateFreeFitParams(p);
            uh = obj.paramStruct;
            obj.simMatrix = diffusion(obj.t,uh.D,uh.gel_radius,uh.finalConc,obj.nmax,obj.rres,obj.rlim,obj.beam_size);
            %obj.simMatrix = obj.simMatrix';
        end
    end
    
end