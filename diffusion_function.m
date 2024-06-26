classdef (Abstract) diffusion_function
    
    
    properties
        %
        %  Things that can't be fit parameters
        %
        
        %time axis properties
        %dt;
        %n_t;
        %n_zp;
        %n_under = 0;
        
        t;
        
        %system properties
        gel_radius;
        beam_size;
        
        dataMatrix;
        
        % fit parameters
        %tolfun relates to the squared residual in the error
        %tolx relates to the displacement in your parameters
        tolfun = 1e-6;%1e-17?
        tolx = 1e-6;%1e-10?
        maxfun = 1e3; %maximum number of function evaluations
        
        %container for fit parameters, these are the ones actually used for fitting
        paramStruct; %structure of all parameters (fixed and free)
        freeParamNames; %names of only free params (cell array)
        p0; %starting point (vector)
        %weightMatrix; %optional for _w functions
        
        %
        % Things that might be fit parameters
        %

        %
        %  Output
        %
        simMatrix;
        spec; %a spectrum just for testing
        pfit; %ending point (vector)
        err; %resulting error
        fitResult; %output structure for the fit result
        
        % from bounded version
        lb; %lower bounds
        ub; %upper bounds
        useParallel=false;
        %nboot;
        %pboot;
    end
    
    methods (Abstract)
        calcCurve;
        fitFunction;
    end
    methods
        function obj = diffusion_function(options)
            if nargin>0
                % assign properties from the input options struct
                props = properties(obj);
                for ii = 1:length(props)
                    if isfield(options,props{ii})
                        obj.(props{ii}) = options.(props{ii});
                    end
                end
                
                % get parameters initialized
                obj = obj.makeParamStruct;
                obj.freeParamNames = obj.freeFitParamNames;
                obj.p0 = obj.freeFitParamInitialValues;
                obj.lb = obj.freeFitParamLowerBounds;
                obj.ub = obj.freeFitParamUpperBounds;
                             
            end
        end
        
        function obj = updateFreeFitParams(obj,p)
            for ii = 1:length(obj.freeParamNames)
                obj.paramStruct.(obj.freeParamNames{ii}) = p(ii);
            end
        end
        
        function obj = makeParamStruct(obj)
            n = obj.allFitParamNames;
            v = num2cell(obj.allFitParamInitialValues);
            obj.paramStruct = cell2struct(v,n,2);
        end
        
        function [out] = freeFitParamInitialValues(obj)
            props = properties(obj);
            out = [];
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    if obj.(props{ii}).isFree
                        out = [out,obj.(props{ii}).value];
                    end
                end
            end
        end
        function [out] = freeFitParamValues(obj)
            props = properties(obj);
            out = [];
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    if obj.(props{ii}).isFree
                        out = [out,obj.(props{ii}).value];
                    end
                end
            end
        end
        function [out] = freeFitParamNames(obj)
            props = properties(obj);
            out = {};
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    for jj = 1:length(obj.(props{ii}))
                        if obj.(props{ii})(jj).isFree
                            n = obj.(props{ii})(jj).name;
                            out = [out,n]; %square brackets make it so cells don't get nested
                        end
                        %what about out = [out n(obj.(props{ii}).isFree]; if
                        %isFree is made to be an array
                    end
                end
            end
        end
        function [out] = allFitParamInitialValues(obj)
            props = properties(obj);
            out = [];
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    
                    out = [out,obj.(props{ii}).value];
                    
                end
            end
        end
        function [out] = allFitParamNames(obj)
            props = properties(obj);
            out = {};
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    n = [obj.(props{ii}).name];
                    out = [out,n]; %the square brackets make it so the cells don't get nested
                end
            end
        end
        
        function prettyPrintFreeParamValues(obj,p)
            %print the parameter vector with associated names
            for ii = 1:length(p)
                fprintf(1,'%20s\t%12f\n',obj.freeParamNames{ii},p(ii));
            end
        end
        % methods from bounded version
        function [out] = freeFitParamLowerBounds(obj)
            props = properties(obj);
            out = [];
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    for jj = 1:length(obj.(props{ii}))
                        if obj.(props{ii})(jj).isFree
                            out = [out,obj.(props{ii})(jj).lb];
                        end
                    end
                end
            end
        end
        function [out] = freeFitParamUpperBounds(obj)
            props = properties(obj);
            out = [];
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    for jj = 1:length(obj.(props{ii}))
                        if obj.(props{ii})(jj).isFree
                            out = [out,obj.(props{ii})(jj).ub];
                        end
                    end
                end
            end
        end
        function obj = globalFit(obj)%dataMatrix,w1,w3,p_array,gfstruct,lb,ub)
            
            % set fit parameters (I believe based on version of matlab)
            if exist('optimoptions','file'),
                opt = optimoptions('fmincon','Algorithm','active-set','Display','iter','TolFun',obj.tolfun,'TolX',obj.tolx,'MaxFunEvals',obj.maxfun);
            elseif exist('optimset','file'),
                %    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
                opt = optimset('Algorithm','active-set','Display','iter','TolFun',obj.tolfun,'TolX',obj.tolx,'MaxFunEvals',obj.maxfun);
            else
                warning('Could not set parameters! Look for the right options functon for your installation.');
            end
            
            opt = optimoptions(opt,'UseParallel',obj.useParallel);
            
            tic
            
            %fmincon has required parameters of error function and initial guess.
            %Documentation has a bunch of additional parameters, most of which we don't
            %understand, but the syntax for not using them is to leave them as blanks.
            [pfit,err] = fmincon(@(p)obj.err_fun(p),obj.p0,[],[],[],[],obj.lb,obj.ub,[],opt);
            
            toc
            
            obj.pfit = pfit;
            obj.err = err;
            
            for ii = 1:length(pfit)
                fprintf(1,'%20s\t%12f\n',obj.freeParamNames{ii},pfit(ii));
            end
        end
        function obj = globalFitBootstrap(obj)%dataMatrix,w1,w3,p_array,gfstruct,lb,ub)
            
            % set fit parameters (I believe based on version of matlab)
            if exist('optimoptions','file'),
                opt = optimoptions('fmincon','Algorithm','active-set','Display','iter','TolFun',obj.tolfun,'TolX',obj.tolx,'MaxFunEvals',obj.maxfun);
            elseif exist('optimset','file'),
                %    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
                opt = optimset('Algorithm','active-set','Display','iter','TolFun',obj.tolfun,'TolX',obj.tolx,'MaxFunEvals',obj.maxfun);
            else
                warning('Could not set parameters! Look for the right options functon for your installation.');
            end
            
            opt = optimoptions(opt,'UseParallel',obj.useParallel);
            
            obj.pboot = zeros(obj.nboot,length(obj.p0));
            npoints = numel(obj.dataMatrix);
            
            
            % go!
            bootStart = tic;
            for ii = 1:obj.nboot;
                ind = randi(npoints,1,npoints);
                %[pboot(ii,:),err] = fmincon(@(p)obj.err_fun(p),obj.p0,[],[],[],[],obj.lb,obj.ub,[],opt);
                [obj.pboot(ii,:)] = fmincon(@(p,ind)obj.err_fun_boot(p,ind),obj.p0,[],[],[],[],obj.lb,obj.ub,[],opt,ind);
            end
            bootEnd = toc(bootStart);
            
            fprintf(1,'Elapsed time: %i seconds\n',bootEnd);
            % the 2 sigma limits should be 95% intervals
            for ii = 1:size(obj.pboot,2)
                fprintf(1,'%20s\t%12f\tpm\t%12.3f\n',...
                    obj.freeParamNames{ii},mean(obj.pboot(:,ii)),2*std(obj.pboot(:,ii)));
            end
            
        end
        function obj = globalFitW(obj)%dataMatrix,w1,w3,p_array,gfstruct,lb,ub)
            
            % set fit parameters (I believe based on version of matlab)
            if exist('optimoptions','file'),
                opt = optimoptions('fmincon','Algorithm','active-set','Display','iter','TolFun',obj.tolfun,'TolX',obj.tolx,'MaxFunEvals',obj.maxfun);
            elseif exist('optimset','file'),
                %    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
                opt = optimset('Algorithm','active-set','Display','iter','TolFun',obj.tolfun,'TolX',obj.tolx,'MaxFunEvals',obj.maxfun);
            else
                warning('Could not set parameters! Look for the right options functon for your installation.');
            end
            
            opt = optimoptions(opt,'UseParallel',obj.useParallel);
            
            tic
            
            %fmincon has required parameters of error function and initial guess.
            %Documentation has a bunch of additional parameters, most of which we don't
            %understand, but the syntax for not using them is to leave them as blanks.
            [pfit,err] = fmincon(@(p)obj.err_fun_w(p),obj.p0,[],[],[],[],obj.lb,obj.ub,[],opt);
            
            toc
            
            obj.pfit = pfit;
            obj.err = err;
            
            for ii = 1:length(pfit)
                fprintf(1,'%20s\t%12f\n',obj.freeParamNames{ii},pfit(ii));
            end
        end
        
        function obj = globalFitBootstrapW(obj)%dataMatrix,w1,w3,p_array,gfstruct,lb,ub)
            
            % set fit parameters (I believe based on version of matlab)
            if exist('optimoptions','file'),
                opt = optimoptions('fmincon','Algorithm','active-set','Display','iter','TolFun',obj.tolfun,'TolX',obj.tolx,'MaxFunEvals',obj.maxfun);
            elseif exist('optimset','file'),
                %    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
                opt = optimset('Algorithm','active-set','Display','iter','TolFun',obj.tolfun,'TolX',obj.tolx,'MaxFunEvals',obj.maxfun);
            else
                warning('Could not set parameters! Look for the right options functon for your installation.');
            end
            
            opt = optimoptions(opt,'UseParallel',obj.useParallel);
            
            obj.pboot = zeros(obj.nboot,length(obj.p0));
            npoints = numel(obj.dataMatrix);
            
            
            % go!
            bootStart = tic;
            for ii = 1:obj.nboot;
                ind = randi(npoints,1,npoints);
                %[pboot(ii,:),err] = fmincon(@(p)obj.err_fun(p),obj.p0,[],[],[],[],obj.lb,obj.ub,[],opt);
                [obj.pboot(ii,:)] = fmincon(@(p,ind)obj.err_fun_boot_w(p,ind),obj.p0,[],[],[],[],obj.lb,obj.ub,[],opt,ind);
            end
            bootEnd = toc(bootStart);
            
            fprintf(1,'Elapsed time: %i seconds\n',bootEnd);
            % the 2 sigma limits should be 95% intervals
            for ii = 1:size(obj.pboot,2)
                fprintf(1,'%20s\t%12f\tpm\t%12.3f\n',...
                    obj.freeParamNames{ii},mean(obj.pboot(:,ii)),2*std(obj.pboot(:,ii)));
            end
            
        end
        % deciding whether or not to keep these methods
        function chi2 = err_fun(obj,p)
            
            obj = obj.calcCurve(p);
            chi2 = sum(sum(sum((obj.dataMatrix-obj.simMatrix).^2)));
            
        end
        function chi2 = err_fun_boot(obj,p,ind)
            
            obj = obj.calcCurve(p);
            chi2 = sum(sum(sum((obj.dataMatrix(ind)-obj.simMatrix(ind)).^2)));
            
        end
       function chi2 = err_fun_w(obj,p)
            
            obj = obj.calcCurve(p);
            chi2 = sum(sum(sum(((obj.dataMatrix-obj.simMatrix).^2).*obj.weightMatrix)));
            
        end
        function chi2 = err_fun_boot_w(obj,p,ind)
            
            obj = obj.calcCurve(p);
            chi2 = sum(sum(sum(((obj.dataMatrix(ind)-obj.simMatrix(ind)).^2).*obj.weightMatrix)));
            
        end
        function out = residuals(obj,p)
            if nargin>1
                 obj = obj.calcCurve(p);
            end
            out = obj.dataMatrix-obj.simMatrix;
        end
        function out = residuals_w(obj,p)
            if nargin>1
                obj = obj.calcCurve(p);
            end
            out = (obj.dataMatrix-obj.simMatrix).*sqrt(obj.weightMatrix);
        end
        function out = residualsSq(obj,p)
            if nargin>1
                obj = obj.calcCurve(p);
            end
            out = (obj.dataMatrix-obj.simMatrix).^2;
        end
        function out = residualsSq_w(obj,p)
            if nargin>1
                obj = obj.calcCurve(p);
            end
            out = (obj.dataMatrix-obj.simMatrix).^2.*obj.weightMatrix;
        end
        
    end
end
