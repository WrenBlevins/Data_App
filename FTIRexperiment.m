classdef FTIRexperiment
    properties
        data double = 0
        freqAxis (:,1) double = 0
        volume (:,1) double = 1 % must be in microliters
        pathLength (:,1) double = 12 % must be in micrometers
        radius (:,1) double % not part of constructor method
        timeInterval (:,1) double = 0 % must be in seconds
        gasFactor (1,1) double
        sample string = ""
        diffusionCoeff double = 0
        finalSpectrum double = 0 % to get the concentration at saturation
        finalConc double = 0 % concentration at saturation
        fittedSpectra struct
    end
    methods
        function f = FTIRexperiment(data,freqAxis,volume,pathLength,timeInterval,sample)
            if nargin > 0
                f.data = data;
                f.freqAxis = freqAxis;
                f.volume = volume;
                f.pathLength = pathLength;
                f.sample = sample;
                f.timeInterval = timeInterval;
                f.radius = sqrt(f.volume*1e9/(pi*f.pathLength)); % returns radius in micrometers
            end
        end
        function plts = plotSpectra(f,specNum,gasFactor)
            % plots the spectra for data f for specNum indices, w optional
            % gasFactor
            %syntax: plotSpectra(f,specNum,gasFactor)
            if numel(specNum) == 1
                if nargin == 3
                    %find the indicies for the amount of spectra desired
                    spectraIndicies = zeros(1,specNum);
                    interval = floor(size(f.data,2)/specNum);
                    for ii = 0:specNum-1
                        spectraIndicies(ii+1) = 1+(ii*interval);
                    end
                    if spectraIndicies(end) ~= size(f.data,2)
                        spectraIndicies = [spectraIndicies size(f.data,2)];
                    end
                    %baseline correction
                    baseline = f.data - f.data(28375,:);
                    %gas line correction
                    %cd("Users/matthewliberatore/Library/CloudStorage/OneDrive-UniversityofPittsburgh/data/ftir_data/Matt/Gas Lines Ref")
                    wha = load("CO2_gas_lines.mat",'gasLines');
                    fixed = baseline-gasFactor.*wha.gasLines;
                elseif nargin == 2
                    %find the indicies for the amount of spectra desired
                    spectraIndicies = zeros(1,specNum);
                    interval = floor(size(f.data,2)/specNum);
                    for ii = 0:specNum-1
                        spectraIndicies(ii+1) = 1+(ii*interval);
                    end
                    if spectraIndicies(end) ~= size(f.data,2)
                        spectraIndicies = [spectraIndicies size(f.data,2)];
                    end
                    %baseline correction
                    baseline = f.data - f.data(28375,:);
                    %gas line correction
                    %cd("Users/matthewliberatore/Library/CloudStorage/OneDrive-UniversityofPittsburgh/data/ftir_data/Matt/Gas Lines Ref")
                    wha = load("CO2_gas_lines.mat",'gasLines');
                    fixed = baseline-f.gasFactor.*wha.gasLines;
                end
            elseif numel(specNum) > 1
                if nargin == 3
                    spectraIndicies = specNum;
                    %baseline correction
                    baseline = f.data - f.data(28375,:);
                    %gas line correction
                    %cd("Users/matthewliberatore/Library/CloudStorage/OneDrive-UniversityofPittsburgh/data/ftir_data/Matt/Gas Lines Ref")
                    wha = load("CO2_gas_lines.mat",'gasLines');
                    fixed = baseline-gasFactor.*wha.gasLines;
                elseif nargin == 2
                    %find the indicies for the amount of spectra desired
                    spectraIndicies = specNum;
                    %baseline correction
                    baseline = f.data - f.data(28375,:);
                    %gas line correction
                    %cd("Users/matthewliberatore/Library/CloudStorage/OneDrive-UniversityofPittsburgh/data/ftir_data/Matt/Gas Lines Ref")
                    wha = load("CO2_gas_lines.mat",'gasLines');
                    fixed = baseline-f.gasFactor.*wha.gasLines;
                end
            end
            plts = plot(f.freqAxis(:,1),fixed(:,spectraIndicies));
            
            hold on
            xlim([2290 2390])
            xlabel("Frequency (cm^{-1})")
            ylabel("Absorbance (AU)")
            legend(string(spectraIndicies))
        end
        function conc = concOverTime(f)
            %baseline correction
            baseline = f.data - f.data(28375,:);
            %gas line correction
            %cd("Users/matthewliberatore/Library/CloudStorage/OneDrive-UniversityofPittsburgh/data/ftir_data/Matt/Gas Lines Ref")
            wha = load("CO2_gas_lines.mat",'gasLines');
            fixed = baseline-f.gasFactor.*wha.gasLines;
            
            CO2band = fixed(f.freqAxis(:,1) > 2290 & f.freqAxis(:,1) < 2390,:);
            conc = max(CO2band)./(1000*f.pathLength*1e-4);
            conc = conc-conc(1); % DOES THIS NEED TO BE HERE?
        end
        function axis = timeAxis(f)
            % generates time axis for data in f depending on time interval
            % specified
            %syntax: timeAxis(f)
            axis = (0:(size(f.data,2)-1)).*10;
        end
        function plts = plotConcOverTime(f)
            %converts data in f to vector of concentration values of CO2
            %for each spectrum, generates time axis, returns a plot of
            %concentration vs time
            %syntax: plotConcOverTime(f)
            plts = plot(timeAxis(f),concOverTime(f));
            hold on
            xlabel("Time (s)")
            ylabel("Concentration (M)")
            hold off
        end
        function f = getFinalConc(f)
            %baseline correction
            baseline = f.finalSpectrum - f.finalSpectrum(28375);
            %gas line correction
            %cd("Users/matthewliberatore/Library/CloudStorage/OneDrive-UniversityofPittsburgh/data/ftir_data/Matt/Gas Lines Ref")
            wha = load("CO2_gas_lines.mat",'gasLines');
            fixed = baseline-f.gasFactor.*wha.gasLines;
            
            CO2band = fixed(f.freqAxis(:,1) > 2290 & f.freqAxis(:,1) < 2390,:);
            f.finalConc = max(CO2band)./(1000*f.pathLength*1e-4);
        end
        function f = gasLineFit(f,center,wg,wl,a1,a2,a3,c0,c1)
            
            n_spectra = size(f.data,2); % number of columns
            
            %initial guess from inputs do this before calling function
            sp = [center,wg,wl,a1,a2,a3,c0,c1];
            %upper and lower bounds
            lb = [2300, 0.5, 0.5,   0, 0.0,   0, -10, -1];
            ub = [2400, 4,   4,   100, 0.2, inf,  10,  1];
            
            opts = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',lb,'Upper',ub,'StartPoint',sp,...
                'Display','Iter');
            ft = fittype(@(center,w_g,w_l,a1,a2,a3,c0,c1,w) co2GasLineFitFunction(w,center,w_g,w_l,a1,a2,a3,c0,c1),...
                'independent',{'w'},'dependent','absorbance',...
                'coefficients',{'center','w_g','w_l','a1','a2','a3','c0','c1'},...
                'options',opts);
            
            %clear out
            out(n_spectra) = struct('x',[],'ydata',[],'yfit',[],'res',[],...
                'fobj',[],'G',[],'O',[]);
            
            % start a timer
            tic
            
            % set the fit range
            range1 = [2290 2390];
            
            % fit each spectrum
            for ii = 1:n_spectra
                
                freq = flip(f.freqAxis);
                s = flip(f.data(:,ii));
                
                % update the fitting region (x and y)
                ind1 = find(freq>=range1(1) & freq<range1(2));
                x = freq(ind1);
                ydata = s(ind1);
                
                % do the fit
                [fobj, G, O] = fit(x,ydata,ft);
                
                
                % get fit result for plotting
                yfit = fobj(x);
                
                % pack up the data and results
                out(ii).x = x;
                out(ii).ydata = ydata;
                out(ii).yfit = yfit;
                out(ii).res = ydata - yfit;
                out(ii).fobj = fobj;
                out(ii).G = G;
                out(ii).O = O;
                
            end
            
            % stop the timer
            toc
            
            % check results
            for ii = 1:n_spectra
                if out(ii).O.exitflag < 1
                    warning('Spectrum %i did not converge!!! Results might not be trustworthy.',ii);
                end
            end
            f.fittedSpectra = out;
        end
        function concs = fittedConcOverTime(f)
            n_spectra = numel(f.fittedSpectra);
            OD = zeros(1,n_spectra);
            for ii = 1:n_spectra
                fobj = f.fittedSpectra(ii).fobj; % copy of the fit
                fobj.a3 = 0; % set gas lines to zero
                fobj.c0 = 0; % set the offset to zero
                fobj.c1 = 0; % set the slope to zero
                
                % call the changed fit (we have to input a vector for it to treat it as an
                % input value to evaluate the fit at)
                temp = fobj([fobj.center fobj.center]);
                OD(ii) = temp(1); %take the first value only (they are duplicates)
            end
            concs = OD./(1000*f.pathLength*1e-4);
            
        end
    end
end