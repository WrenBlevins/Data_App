function out = diffusion_moving_beam(t,D,A,C,nmax,sigma,dx,dy,varargin)
% out = diffusion_moving_beam(t,D,A,C,nmax,sigma,dx,dy,varargin)
% DIFFUSION FITTING FUNCTION. varargin currently only has two options: rlim and t0.
% Enter as name/value pairs.
% The default values reflect the case where rlim >> sigma (Inf) and t0 = 0.

j0 = besselzero(0,nmax,1); % get nmax zeros of 0th order bessel function of first kind

c0 = zeros(1,nmax); % constants defined by boundary condition
for ii = 1:nmax
    c0(ii) = (C*2)/(j0(ii)*besselj(1,j0(ii)));
end

x = linspace(-A,A,50);
y = linspace(-A,A,50);

if mod(numel(varargin),2) ~= 0
    error("Watch for unpaired name/value pairs!")
end
while numel(varargin) >= 2
    var = varargin{1};
    val = varargin{2};
    switch var
        case "rlim"
            rlim = val;
        case "t0"
            t0 = val;
        otherwise
            error("Invalid name/value pair")
    end
    varargin = varargin(3:end);
end
if ~exist('rlim','var')
   rlim = 1000*sigma; 
end
if ~exist('t0','var')
    t0 = 0;
end

if isscalar(t)
%     u = zeros(1,numel(r));
%     for ii = 1:nmax
%         u = u + c0(ii).*besselj(0,j0(ii)/A.*r).*exp(-(j0(ii)/A)^2*D.*t);
%     end
%     u = (C-u).*exp(-r.^2/(2*sigma^2)).*r*...
%         (-1/(sigma^2*(exp(-rlim^2/(2*sigma^2)) - 1)));
%     out = trapz(r,u);
      error("Time axis parameter t must be a vector, not a scalar.")
else
    [X,Y] = meshgrid(x,y);
    mask = ones(size(X));
    mask(X.^2+Y.^2 > rlim.^2) = 0;
    norm = trapz(y,trapz(x,exp(-((X-dx).^2+(Y-dy).^2)/(2*sigma^2)).*mask,1));
    clear X Y
    
    [X,Y,T] = ndgrid(x,y,t+t0);
    u = zeros(numel(x),numel(y),numel(t));
    for ii = 1:nmax
        u = u + c0(ii).*besselj(0,j0(ii)/A.*sqrt(X.^2+Y.^2)).*exp(-(j0(ii)/A)^2*D.*T);
    end
    u = (C-u).*exp(-((X-dx).^2+(Y-dy).^2)/(2*sigma^2))/norm;
    u = u.*mask;
    out = squeeze(trapz(y,trapz(x,u,1),2));
end
end