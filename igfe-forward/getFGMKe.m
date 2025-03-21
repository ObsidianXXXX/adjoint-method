function [keVec,ke] = getFGMKe(Dgauss,Bgauss,detJgauss,w)
%%%===========================Copyright==================================%%%
	%%%   Version Mar. 2024
	%%%
	%%%   Gengxuan Zhu <12324117@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This function is to calculate the stiffness matrix of each element  
	%%% at each Gauss Point.
	%%%======================================================================%%%
	%%%======================================================================%%%
    for i = 1:length(w)
		Bgauss_i = Bgauss(:,i); detJgauss_i = detJgauss(:,i);
        DGauss_i = Dgauss(:,i);
		BDBJgauss_i = cellfun(@(B,detJ,stiff)B'*stiff*B*detJ,Bgauss_i,...
			detJgauss_i,DGauss_i,'UniformOutput',0);
	    wi = w(i);
		if i==1 
			ke = cellfun(@(BDBJ)wi*BDBJ,BDBJgauss_i,'UniformOutput',0);
        else
			ke = cellfun(@(k,BDBJ)k+wi*BDBJ,ke,BDBJgauss_i,'UniformOutput',0);
        end
    end
    keVec =cellfun(@(Ke)Ke(tril(true(size(Ke)))),ke,'UniformOutput',0);
end