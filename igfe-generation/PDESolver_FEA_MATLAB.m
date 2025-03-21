function state = PDESolver_FEA_MATLAB(feInfo,Eele,bcInfo);
%%%===========================Copyright==================================%%%
%%%   Xueyan Hu <huxueyan@zju.edu.cn>
%%%   Institute of Applied Mechanics,Zhejiang University
%%%
%%%===========================Description================================%%%
%%% A matlab FEA solver
%%% Lastest version May. 2021
%%% ----INPUT
%%% meshInfo: <structure> mesh information
%%% materialInfo: <structure> material information
%%% BCInfo:<structure> boundary condition information
%%% preMatrixInfo:<structure> boundary condition information
%%%
%%% ----OUTPUT
%%% FEAResult: <structure>
%%%
%%%======================================================================%%%
	% BC info
	fixdof = bcInfo.fixdof; Fall = bcInfo.F;
	if iscell(Fall)
		for i = 1:length(Fall)
			F = Fall{i};
			FEAResult = solve_single(feInfo,Eele,fixdof,F);
			if length(Fall)==1
				state = FEAResult;
			else
				state{i} = FEAResult;
			end
		end
	else
		state = solve_single(feInfo,Eele,fixdof,Fall);
	end
	
		

%% -- Solve for one single case;
function FEAResult = solve_single(feInfo,Eele,fixdof,F)
	meshInfo = feInfo.meshInfo; materialInfo = feInfo.materialInfo;
	preMatrixInfo = feInfo.preMatrixInfo; 
	% Mesh Info
	Ivar = meshInfo.Ivar; nEl = meshInfo.nEl;
	eleDofID = meshInfo.eleDofID;
	% K info
	Ke0Vec = preMatrixInfo.KeVec;
	% K = sparse(2*nodeSize, 2*nodeSize);
	alldofs = length(F);
	U = zeros(alldofs,1);
	EeleCell = mat2cell(Eele,ones(size(Eele)),1);
	penK = cell2mat(cellfun(@(kVec,E)E*kVec,Ke0Vec,EeleCell,'UniformOutput',0));
	K = fsparse(Ivar(:,1),Ivar(:,2),penK,[alldofs,alldofs]);
	freedofs    = setdiff([1:alldofs],fixdof);
	U(freedofs,:) = decomposition(K(freedofs,freedofs),'chol','lower') \ F(freedofs,:);
    K_tot = K + K' - diag(diag(K));
    FEAResult.RF = K_tot * U;
    if (isfield(feInfo,'outputRequest')==0)
		FEAResult.U = U;
	else
		OutputRequest = feInfo.outputRequest;
		FEAResult.U = U;
		%%-- History output
		if isfield(OutputRequest,'History');
			History = OutputRequest.History;
			for i = 1:length(History)
				if strcmp(History{i},'Comp')
					FEAResult.History.Comp = F'*U;
				elseif strcmp(History{i},'GloK')
					FEAResult.History.GloK = sparse(K+K'-diag(K));
				end
			end
		end

		%%-- Field output at Element
		if isfield(OutputRequest,'FieldAtEle');
			FieldAtEle = OutputRequest.FieldAtEle;
			Ue = U(eleDofID);
			[nel,nDof] = size(Ue);
			ind = find(tril(ones(nDof)));
			Ke0 = cellfun(@(keVec)fullKe(keVec,ind,nDof),Ke0Vec,'UniformOutput',0);
			UeCell = mat2cell(Ue,ones(nel,1),nDof);
			for i = 1:length(FieldAtEle)
				if strcmp(FieldAtEle{i},'EleComp')
					Elecomp = Eele.*cell2mat(cellfun(@(Ue,Ke)Ue*Ke*Ue',UeCell,Ke0,'UniformOutput',0));
					FEAResult.FieldAtEle.EleComp = Elecomp;
				elseif strcmp(FieldAtEle{i},'EleCompSolid')
					ElecompSolid = cell2mat(cellfun(@(Ue,Ke)Ue*Ke*Ue',UeCell,Ke0,'UniformOutput',0));
					FEAResult.FieldAtEle.EleCompSolid = ElecompSolid;
				end
			end
		end

		%%-- Field output at Node
		if isfield(OutputRequest,'FieldAtNode');
			FieldAtNode = OutputRequest.FieldAtNode;
			for i = 1:length(FieldAtNode)
				requestI = FieldAtNode{i};
				pos = requestI{1};
				if strcmp(pos,'At_Gauss')
					B_At_Gauss = preMatrixInfo.B_At_Gauss;
				elseif strcmp(pos,'At_EleCen')
					B_At_EleCen = preMatrixInfo.B_At_EleCen;
				end
			end
			for i = 1:length(FieldAtNode)
				requestI = FieldAtNode{i};
				pos = requestI{1};
				Var = requestI{2};
				calUtoE = sum(find(strcmp(Var,'UtoE')));
				calE = sum(find(strcmp(Var,'Strain')));
				calS = sum(find(strcmp(Var,'Stress')));
				if (calUtoE  ~= 0)
					if strcmp(pos,'At_Gauss')
						FEAResult.FieldAtNode.UtoE_At_Gauss = B_At_Gauss;
					elseif strcmp(pos,'At_EleCen')
						FEAResult.FieldAtNode.UtoE_At_EleCen = B_At_EleCen;
					end
				end
				if (calS ~= 0) || (calE ~= 0)
					if strcmp(pos,'At_Gauss')
						Strain = calStrain(B_At_Gauss,U,eleDofID,pos,nEl);
					elseif strcmp(pos,'At_EleCen')
						Strain = calStrain(B_At_EleCen,U,eleDofID,pos,nEl);
					end
					if (calE ~= 0)
						if strcmp(pos,'At_Gauss')
							FEAResult.FieldAtNode.Strain_At_Gauss = Strain;
						elseif strcmp(pos,'At_EleCen')
							FEAResult.FieldAtNode.Strain_At_EleCen = Strain;
						end
					end
				end
				if (calS ~= 0)
					% Stress0 = calStress(Strain,ones(size(Eele)),D0,nEl);
					if strcmp(pos,'At_Gauss')
						D0 = preMatrixInfo.D_At_Gauss;
                        Stress = calStress(Strain,Eele,D0,nEl);
                        FEAResult.FieldAtNode.Stress_At_Gauss = Stress;
						% FEAResult.FieldAtNode.Stress0_At_Gauss = Stress0;
					elseif strcmp(pos,'At_EleCen')
						D0 = preMatrixInfo.D_At_EleCen;
                        Stress = calStress(Strain,Eele,D0,nEl);
                        FEAResult.FieldAtNode.Stress_At_EleCen = Stress;
						% FEAResult.FieldAtNode.Stress0_At_EleCen = Stress0;
					end
				end
			end
		end
	end

%% -- calculation of Strain
function Strain = calStrain(B,U,eleDofID,pos,nEl)
	Ue = U(eleDofID);
	UeCell = mat2cell(Ue,ones(nEl,1));
	% size(UeCell{1})
	% size(B{1})
	if strcmp(pos,'At_EleCen')
		Strain = cellfun(@(x,y)x*y',B,UeCell,'UniformOutput',0);
	elseif strcmp(pos,'At_Gauss')
		for j = 1:4
			Strain(:,j) =  cellfun(@(x,y)x*y',B(:,j),UeCell,'UniformOutput',0);
		end
	end
%% -- calculation of stress
function Stress = calStress(Strain,Eele,D0,nEl)
	% EeleCell = mat2cell(Eele,ones(nEl,1));
	% DCell = cellfun(@(x,D0)x*D0, EeleCell,D0,'UniformOutput',0);
	Stress = cell(size(Strain));
	for j =  1:length(Strain(1,:))
		Stress(:,j) = cellfun(@(x,y)x*y,D0(:,j),Strain(:,j),'UniformOutput',0);
	end

function y = fullKe(x,ind,nDof)
	yTri = zeros(nDof);
	yTri(ind) = x;
	y = yTri+yTri'-diag(diag(yTri));
