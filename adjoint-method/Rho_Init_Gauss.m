function init_rho = Rho_Init_Gauss(meshInfo,upl,lowl,sigma)
%%%===========================Copyright==================================%%%
	%%%   Version Dec. 2024
	%%%
	%%%   Gengxuan Zhu <zhugx@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a  function to initialize modulus vector in Gauss
    %%% distribution.
	%%% ----INPUT
	%%% meshInfo: <structure> mesh information gennerated form
	%%%				fun{Init_PDEslover_Mesh_UniforRec}
	%%% upl: Upper limit of modulus values at each node.
	%%% lowl: Lower limit of modulus values at each node.
	%%% sigma: Mean value of Gauss distribution.
	%%%
	%%% ----OUTPUT
	%%% init_Rho: Initialized density vector.
	%%%======================================================================%%%
    if nargin < 2
        sigma = 0.5;
        upl = 2;
        lowl = 1;
    elseif nargin < 4
        sigma = 0.5;
    end
    nEl = meshInfo.nNod;
    ave = (upl+lowl)/2;
    init_rho = sigma*randn(nEl,1)+ave;
    while any(init_rho < lowl | init_rho > upl)
        idx = init_rho < lowl | init_rho > upl;
        init_rho(idx) = sigma * randn(sum(idx), 1) + ave;
    end
end