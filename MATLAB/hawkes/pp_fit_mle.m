function [vPar_est, fLle, exitflag, mHessian] = pp_fit_mle(negLLfun,PP,OptArgs)
%
%   vCat  - the catalogue of events
%   Hawkes - input model structure: fT0, fTs, fTe, fT1, fM0, ...
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 31 October 2024
%   ...
%   version: 1.0.0, 2 November 2024
%
    arguments
        negLLfun                              % negative log-likelihood function 
        PP struct                             % structure that defines the poit-process to fit
        OptArgs.OptimMethod char = 'fmincon'  % 'gs'; 'fminsearch';  
        OptArgs.NonLinCon struct = []         % 
        OptArgs.Display char = 'final'        % 'iter';    % fmincon level of display
    end
    
    if ~isempty(OptArgs.NonLinCon)
        %beta = 2.3;
        %brmax = 1.0;
        brcon = @(x) br_nonlcon(OptArgs.NonLinCon.beta,x,OptArgs.NonLinCon.brmax);
    end

    mHessian = zeros(length(PP.mIP.vInitParams));
    if strcmp(OptArgs.OptimMethod,'fmincon')     % problem setup for fmincon uses optimoptions() from Global search toolbox
        %options = optimoptions('fmincon','Algorithm','sqp','StepTolerance',1e-8,'Display',OptArgs.Display);
        %options = optimoptions('fmincon','Algorithm','sqp','StepTolerance',1e-7,'Display',OptArgs.Display,'Diagnostics','on','DiffMinChange',0.01);
        %options = optimoptions('fmincon','Algorithm','sqp','MaxIterations',1000,'MaxFunctionEvaluations',1500,'StepTolerance',1e-8,'Display',OptArgs.Display);
        options = optimoptions('fmincon','Algorithm','sqp','MaxIterations',1000,'MaxFunctionEvaluations',3500,'StepTolerance',1e-10,'Display',OptArgs.Display);
%         problem = createOptimProblem('fmincon','objective',llfun,'x0',ETAS.mIP.vInitParams,...
%                             'lb',ETAS.mIP.vLowB,'ub',ETAS.mIP.vUppB,'Aeq',ETAS.mIP.Aeq,'beq',ETAS.mIP.beq,'options',options);
        problem.solver    = 'fmincon';
        problem.objective = negLLfun;
        problem.x0        = PP.mIP.vInitParams;
        problem.options   = options;
        problem.lb        = PP.mIP.vLowB;
        problem.ub        = PP.mIP.vUppB;
        problem.Aeq       = PP.mIP.Aeq;
        problem.beq       = PP.mIP.beq;
        if ~isempty(OptArgs.NonLinCon)
            problem.nonlcon = brcon;
        end

        [vPar_est,fLle,exitflag,output,lambda,grad,mHessian] = fmincon(problem);
    elseif strcmp(OptArgs.OptimMethod,'fminsearch')
        %options = optimset('Display',sDisplay,'MaxFunEvals',2000,'TolX',1e-7,'TolFun',1e-7);
        options = optimset('Display',OptArgs.Display,'TolX',1e-7,'TolFun',1e-7);
        %options = optimset('Display',sDisplay,'TolX',1e-5,'TolFun',1e-5,'PlotFcns',@optimplotfval);
        [vPar_est,fLle,exitflag,output] = fminsearch(negLLfun,PP.mIP.vInitParams,options);
    elseif strcmp(OptArgs.OptimMethod,'gs')        % problem setup for global search
        options = optimoptions('fmincon','Algorithm','sqp','StepTolerance',1e-7,'Display',OptArgs.Display);
        problem = createOptimProblem('fmincon','objective',negLLfun,'x0',PP.mIP.vInitParams,...
                            'lb',PP.mIP.vLowB,'ub',PP.mIP.vUppB,'Aeq',PP.mIP.Aeq,'beq',PP.mIP.beq,'options',options);
        gs = GlobalSearch;
        [vPar_est,fLle,exitflag,output] = run(gs,problem);
    elseif strcmp(OptArgs.OptimMethod,'ga')        % problem setup for genetic algorithm
        nvars = length(PP.mIP.vInitParams);
        options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf);
        [vPar_est,fLle,exitflag,output] = ga(negLLfun,nvars,[],[],PP.mIP.Aeq,PP.mIP.beq,PP.mIP.vLowB,PP.mIP.vUppB,[],options);
    elseif strcmp(OptArgs.OptimMethod,'sa')        % problem setup for simulated annealing algorithm
        options = optimoptions('simulannealbnd','FunctionTolerance',1e-6,'PlotFcn', @saplotbestf);
        [vPar_est,fLle,exitflag,output] = simulannealbnd(negLLfun,PP.mIP.vInitParams,PP.mIP.vLowB,PP.mIP.vUppB,options);
    end

    fLle = -fLle; % change the sign to reflect that the maximum of the log-likelihood function is needed
end

