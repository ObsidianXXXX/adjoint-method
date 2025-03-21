function femInfo = FEM_Assemble(meshInfo,materialInfo,BCInfo)
%%%===========================Copyright==================================%%%
	%%%   Version Mar. 2024
	%%%
	%%%   Gengxuan Zhu <12324117@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This is a  function to get some matrix value of the GFE model
	%%%======================================================================%%%
    [Bgauss,w,detJgauss,Dgauss,D_init] = getFGMBD(meshInfo,materialInfo);
	[keVec,Ke] = getFGMKe(Dgauss,Bgauss,detJgauss,w);
    femInfo.keVec = keVec;
    femInfo.Ke = Ke;
    femInfo.Bgauss = Bgauss;
    femInfo.detJgauss = detJgauss;
    femInfo.D0 = D_init;
    femInfo.BCInfo = BCInfo; 
    femInfo.meshInfo = meshInfo;
end