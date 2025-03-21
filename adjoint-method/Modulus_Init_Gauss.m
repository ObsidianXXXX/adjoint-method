function init_Evec = Modulus_Init_Gauss(meshInfo,upl,lowl,sigma)
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
	%%% init_Evec: Initialized modulus vector.
	%%%======================================================================%%%
    if nargin < 2
        sigma = 0.5;
        upl = 4;
        lowl = 2;
    elseif nargin < 4
        sigma = 0.5;
    end
    nEl = meshInfo.nNod;
    ave = (upl+lowl)/2;
    init_Evec = sigma*randn(nEl,1)+ave;
    while any(init_Evec < lowl | init_Evec > upl)
        idx = init_Evec < lowl | init_Evec > upl;
        init_Evec(idx) = sigma * randn(sum(idx), 1) + ave;
    end
end