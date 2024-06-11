classdef u1Mdl < diffusion_function
    properties
        D
        finalConc
        rres
        rlim
    end
    methods
        function obj = u1Mdl(options)
            obj@diffusion_function(options);
        end
        function obj = calcCurve(obj,p)
            obj = obj.updateFreeFitParams(p);
            uh = obj.paramStruct;
            obj.simMatrix = variable_1D_diffusion(obj.t,1,@(x,t,u,dudx) uh.D*u,...
                obj.finalConc,uh.rres,uh.rlim,uh.beam_size);
        end
    end
end