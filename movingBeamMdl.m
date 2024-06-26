classdef movingBeamMdl < diffusion_function
    
    properties
        nmax
        rlim
        finalConc
        dx
        D
        t0
    end
    
    methods
        function obj = movingBeamMdl(options)
            obj@diffusion_function(options);
        end
        function obj = calcCurve(obj,p)
            obj = obj.updateFreeFitParams(p);
            uh = obj.paramStruct;
            obj.simMatrix = diffusion_moving_beam(obj.t,uh.D,uh.gel_radius,...
                uh.finalConc,uh.nmax,uh.beam_size,uh.dx,0,'rlim',uh.rlim,...
                't0',uh.t0);
            %obj.simMatrix = obj.simMatrix';
        end
        function obj = fitFunction(obj)
            
        end
    end
    
end