function out = variable_1D_diffusion(t,m,D,c0,rres,rlim,sigma)

arguments
   t (1,:) double
   m (1,1)
   D function_handle
   c0 (1,1) double
   rres
   rlim
   sigma
end

x = 0:rres:rlim;

test = func2str(D);
if ~startsWith(test,'@(x,t,u,dudx)')
   error("Please type the function handle in the form '@(x,t,u,dudx)'")
   return
end

sol = pdepe(m,@heatcyl,@heatic,@heatbc,x,t);

average = 1; %this is in case I ever want to turn off the radius averaging and see the full surface
if average
    %sigma = 0.997;
    N = trapz(x,exp(-x.^2/(2*sigma^2)).*x,2);
    new = sol.*exp(-x.^2/(2*sigma^2)).*x*1/N;
    out = trapz(x,new,2);
else
    out = sol;
end

    function [c,f,s] = heatcyl(x,t,u,dudx)
        c = 1;
        f = D(x,t,u,dudx)*dudx;
        s = 0;
    end
    function u0 = heatic(x)
        %n = 2.404825557695773;
        %u0 = besselj(0,n*x);
        u0 = 0; % initial distribution of x at t = 0
    end
    function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t)
        pl = 0; %ignored by solver since m=1
        ql = 0; %ignored by solver since m=1
        pr = ur-c0; % ur - time distribution at outer edge
        qr = 0;
    end
end