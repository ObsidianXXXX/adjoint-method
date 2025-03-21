function forwU = Forward_Solver(femInfo)
%%%===========================Copyright==================================%%%
	%%%   Version Dec. 2024
	%%%
	%%%   Gengxuan Zhu <12324117@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a  function to get displacement results under load
    %%% condition.
	%%%======================================================================%%%
    % Assemble Total Stiffness Matrix
    meshInfo = femInfo.meshInfo;
    nDof = meshInfo.nDof;
    Ivar = meshInfo.Ivar;
    keVec = cell2mat(femInfo.keVec);
    K = fsparse(Ivar(:,1),Ivar(:,2),keVec,[nDof,nDof]);
    % Initialize Displacement Vector
    forwU = zeros(nDof,1);
    % FEM Solve
    BCInfo = femInfo.BCInfo;
    F = BCInfo.F; fixdof = BCInfo.fixdof;
    freedof = setdiff(1:nDof,fixdof);
    forwU(freedof,:) = decomposition(K(freedof,freedof),'chol','lower') \ F(freedof,:);
end