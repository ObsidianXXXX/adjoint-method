function [Evec_next,grad_next,nextLambda,U_next,iter] = Wolfe_Line_Search(U,forwU,grad,H,femInfo, ...
    materialInfo,E_min,E_max,alpha0,c1,c2)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
    alpha = alpha0;
    max_iter = 1000;
    success = false;
    % Initial Functiona Value and Gradient
    norm_U = U'*U;
    tarFun = ((U-forwU)'*(U-forwU))/2/norm_U;
    iter = materialInfo.iter;
    iter_Evec = materialInfo.E_hisvec{iter};
    iter = iter+1;
    for i = 1:max_iter
        %% Calculate new Function Value and Gradient
        dk = -H*grad;
        Evec_next = iter_Evec+alpha*dk;
        % Limit Modulus Value
        minus_idx = find(Evec_next<0);
        max_idx = find(Evec_next>E_max);
        Evec_next(minus_idx) = E_min;
        Evec_next(max_idx) = E_max;
        % Update Material Information
        materialInfo = Material_Update(materialInfo,Evec_next,iter);
        % Assemble FEM Information
        femInfo = FEM_Assemble(femInfo.meshInfo,materialInfo,femInfo.BCInfo);
        % Solve the Forward Problem
        U_next = Forward_Solver(femInfo);
        % Functiona Value
        tarFun_next = ((U-U_next)'*(U-U_next))/2/norm_U;
        nextLambda = Lambda_Distance(femInfo,U_next,U);
        grad_next = Get_Grad(femInfo,nextLambda,U_next);
        
        if tarFun_next > tarFun+c1*alpha*(grad'*dk)
            alpha = alpha/2;
        elseif grad_next'*dk < c2*grad'*dk
            alpha = alpha*2;
        else
            success = true;
            break;
        end 
    end
    if ~success
            warning('Wolfe line search did not converge');
    end
end