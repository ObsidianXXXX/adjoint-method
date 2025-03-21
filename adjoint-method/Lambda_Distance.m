function lambda = Lambda_Distance(femInfo,forwU,U)
%%%===========================Copyright==================================%%%
	%%%   Version Dec. 2024
	%%%
	%%%   Gengxuan Zhu <12324117@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a  function to get Lagarange multipliers for distance
    %%% funciton.
	%%%======================================================================%%%
    % Assemble Total Stiffness Matrix
    meshInfo = femInfo.meshInfo;
    nDof = meshInfo.nDof;
    Ivar = meshInfo.Ivar;
    keVec = cell2mat(femInfo.keVec);
    K = fsparse(Ivar(:,1),Ivar(:,2),keVec,[nDof,nDof]);
    % Initialize Displacement Vector
    lambda = zeros(nDof,1);
    % FEM Solve
    norm_U = U'*U;
    delta = (forwU-U)/norm_U;
    BCInfo = femInfo.BCInfo;
    fixdof = BCInfo.fixdof;
    freedof = setdiff(1:nDof,fixdof);
    lambda(freedof,:) = -decomposition(K(freedof,freedof),'chol','lower') \ delta(freedof,:);
end