function materialInfo = Material_Update(materialInfo,iter_Evec,iter)
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
        materialInfo.E_inivec = iter_Evec;
        materialInfo.E_hisvec = {};
    end
    materialInfo.iter = iter;
    materialInfo.E_hisvec{iter} = iter_Evec;
end