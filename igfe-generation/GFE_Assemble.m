function [feInfo] = GFE_Assemble(meshInfo,elementType,materialInfo,...
	outputRequest)
	%%%===========================Copyright==================================%%%
	%%%   Version Mar. 2024
	%%%
	%%%   Gengxuan Zhu <12324117@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%   Xueyan Hu <huxueyan@zju.edu.cn>
	%%%   Qihang Union & Innovation Center (QUIC), Huanjiang Laboratory, Zhuji 311800, China
	%%%
	%%%===========================Description================================%%%
	%%% This is a  function to get some matrix value of the GFE model
	%%% ----INPUT
	%%% meshInfo: <structure> mesh information gennerated form
	%%%				fun{Init_PDEslover_Mesh_UniforRec}
	%%% elementType: Type of the element only CPS4 is avaliable.
	%%% materialInfo: <structure> material information only linear elastic is
	%%%				avaliable
	%%% outRequest: <structure> if stress and strain is request the UtoE map
	%%%				and D0 will be included in the output
	%%%
	%%% ----OUTPUT
	%%% preMatrixInfo: <structure>
	%%% 	Basic Information:
	%%%			preMatrixInfo.Ke0: element stiffness matrix with E=1
	%%%		Extension Information
	%%%			preMatrixInfo.UtoE: matrix to map displacement
	%%%											to strian for each element
	%%%			preMatrixInfo.UtoE: matrix to map strain
	%%%											to stress for each element
	%%%======================================================================%%%
	if isfield(outputRequest,'FieldAtNode')
		cen = 1;
	else
		cen = 0;
	end
	[Bgauss,w,detJgauss,Bcen,Dgauss,Dcen,materialInfo] = ...
				getFGMBD(meshInfo,materialInfo,elementType,cen);
	[ke,keVec,ve] = getFGMKe(meshInfo,Dgauss,Bgauss,detJgauss,w);

	preMatrixInfo.KeVec = keVec;
    preMatrixInfo.Ke = ke;
	preMatrixInfo.ve = cell2mat(ve);

	if isfield(outputRequest,'FieldAtNode')
		position = {};
		Fan = outputRequest.FieldAtNode;
		for i = 1:length(Fan)
			pos = Fan{i}{1};
			if strcmp(pos,'At_Gauss')
				preMatrixInfo.B_At_Gauss = Bgauss;
				preMatrixInfo.D_At_Gauss = Dgauss;
			elseif strcmp(pos,'At_EleCen')
				preMatrixInfo.B_At_EleCen = Bcen;
				preMatrixInfo.D_At_Elecen = Dcen;
			end
		end
    end

	preMatrixInfo.eleCen = PDEInfo_ELeCen(meshInfo);

	meshInfo.elementType = elementType;
	
	feInfo.preMatrixInfo = preMatrixInfo; feInfo.meshInfo = meshInfo;
	feInfo.materialInfo = materialInfo; feInfo.outputRequest = outputRequest;
