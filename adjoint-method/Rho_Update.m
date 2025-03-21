function rhoInfo = Rho_Update(rhoInfo,iter_rho,iter)
%%%===========================Copyright==================================%%%
	%%%   Version Dec. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a  function to update modulus values after iterations.
	%%% ----INPUT
	%%% iter_Evec: Modulus vector after each iteration.
	%%% iter: Iteration numbers.
	%%%
	%%% ----OUTPUT
	%%% materialInfo: Information of material properties.
    %%%     materialInfo.iter: Total Iteration numbers.
    %%%     materialInfo.E_inivec: Initialized modulus vector.
    %%%     materialInfo.E_hisvec: Modulus vector values of all iterations.
	%%%======================================================================%%%
    if nargin < 3
        iter = 1;
        rhoInfo.Rho_inivec = iter_rho;
        rhoInfo.Rho_hisvec = {};
    end
    rhoInfo.iter = iter;
    rhoInfo.Rho_hisvec{iter} = iter_rho;
end