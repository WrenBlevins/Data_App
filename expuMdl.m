classdef expuMdl < diffusion_function
    properties
        D
        b
        f0
        finalConc
        rres
        rlim
    end
    methods
        function obj = expuMdl(options)
            obj@diffusion_function(options);
        end
        function obj = calcCurve(obj,p)
            obj = obj.updateFreeFitParams(p);
            uh = obj.paramStruct;
            obj.simMatrix = variable_1D_diffusion(obj.t,1,...
                @(x,t,u,dudx) uh.D*exp(-u/uh.b)+uh.f0,...
                obj.finalConc,uh.rres,uh.rlim,uh.beam_size);
        end
    end
end