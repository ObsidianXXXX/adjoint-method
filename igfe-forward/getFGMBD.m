function [Bgauss,w,detJgauss,Dgauss,D_int] = getFGMBD(meshInfo,materialInfo)
%%%===========================Copyright==================================%%%
	%%%   Version Mar. 2024
	%%%
	%%%   Gengxuan Zhu <12324117@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This function is to get the FGM element and stiffness matrix of FGM 
	%%% at Gauss Point based on Isoparametric Elements.
	%%%======================================================================%%%
	%%%======================================================================%%%

    [N,dNs,dNt,w] = shapeFunAtGauss();
	elenodsID = meshInfo.eleNodsID;
	nEle = meshInfo.nEl;
	xEle = meshInfo.xEle; yEle = meshInfo.yEle;
    %% Assign Modulus Value for each element
    E_vec = materialInfo.E_hisvec{materialInfo.iter};
    GEle = E_vec(elenodsID);
    nu = materialInfo.nu;

	% Get FGM element
	[nIntPoint,nShapeFun] = size(N);
	Bgauss = []; detJgauss =[];
	Dgauss = [];
	D_int = [1 , nu , 0;
		 nu, 1  , 0;
		 0 , 0  , 0.5*(1-nu)]/(1-nu^2);
	for i = 1:nIntPoint
		dNs_i = dNs(i,:); dNt_i = dNt(i,:);
		%% - J at gauss_i
		dxds = dNs_i*xEle'; dxdt = dNt_i*xEle';
		dyds = dNs_i*yEle'; dydt = dNt_i*yEle'; %size:(1*Nele)
		detJ = dxds.*dydt-dyds.*dxdt;
		dxds = dxds./detJ; dxdt = dxdt./detJ;
		dyds = dyds./detJ; dydt = dydt./detJ;
		%% - grad(u) & grad(v) in s-t coord-
		% dN*a; with a={u1,v2,...un,vn}
		duds = zeros(1,nShapeFun*2); dudt = duds;
		dvds = duds; dvdt = dvds;
		duds(1,1:2:end) = dNs_i; dudt(1,1:2:end) = dNt_i;
		dvds(1,2:2:end) = dNs_i; dvdt(1,2:2:end) = dNt_i;
		%% - grad(u) & grad(v) in x-y coord-
		dudx = (dydt)'*duds-(dyds)'*dudt;
		dudy = -(dxdt)'*duds+(dxds)'*dudt;
		dvdx = (dydt)'*dvds-(dyds)'*dvdt;
		dvdy = -(dxdt)'*dvds+(dxds)'*dvdt;
		%% - get strain
		ex = dudx; ey = dvdy;
		exy = dudy+dvdx;
		%% - Get B at gauss_i
		B(1,:,:) = (ex)';
		B(2,:,:) = (ey)';
		B(3,:,:) = (exy)';
		%% - 3D mat B to cell B
		Bele = mat2cell(B,[3],[8],ones(1,nEle));
		Bele_gauss_i = (Bele(:));
		Bgauss = [Bgauss,Bele_gauss_i];
		detJ_gauss_i = mat2cell(detJ',ones(1,nEle),1);
		detJgauss = [detJgauss,detJ_gauss_i];
		 
		N_i = N(i,:);
		Eele_gauss_i = N_i*GEle';
		D = kron(Eele_gauss_i',D_int);
		Dele = mat2cell(D,3*ones(1,nEle),[3]);
		Dgauss = [Dgauss,Dele];
    end
end