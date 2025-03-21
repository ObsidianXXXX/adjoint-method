function [ke,keVec,ve] = getFGMKe(meshInfo,Dgauss,Bgauss,detJgauss,w)
	nEl = meshInfo.nEl;
	if iscell(Dgauss)
		for i = 1:length(w)
			Bgauss_i = Bgauss(:,i); detJgauss_i = detJgauss(:,i);
            DGauss_i = Dgauss(:,i);
			BDBJgauss_i = cellfun(@(B,detJ,stiff)B'*stiff*B*detJ,Bgauss_i,...
				detJgauss_i,DGauss_i,'UniformOutput',0);
			wi = w(i);
			if i==1
				ke = cellfun(@(BDBJ)wi*BDBJ,BDBJgauss_i,'UniformOutput',0);
				ve = cellfun(@(detJ)wi*detJ,detJgauss_i,'UniformOutput',0);
			else
				ke = cellfun(@(k,BDBJ)k+wi*BDBJ,ke,BDBJgauss_i,'UniformOutput',0);
				ve = cellfun(@(v,detJ)v+wi*detJ,ve,detJgauss_i,'UniformOutput',0);
			end
		end
	else
		for i = 1:length(w)
			Bgauss_i = Bgauss(:,i); detJgauss_i = detJgauss(:,i); 
			BDBJgauss_i = cellfun(@(B,detJ)B'*Dgauss*B*detJ,Bgauss_i,...
				detJgauss_i,'UniformOutput',0);
			wi = w(i);
			if i==1
				ke = cellfun(@(BDBJ)wi*BDBJ,BDBJgauss_i,'UniformOutput',0);
				ve = cellfun(@(detJ)wi*detJ,detJgauss_i,'UniformOutput',0);
			else
				ke = cellfun(@(k,BDBJ)k+wi*BDBJ,ke,BDBJgauss_i,'UniformOutput',0);
				ve = cellfun(@(v,detJ)v+wi*detJ,ve,detJgauss_i,'UniformOutput',0);
			end
		end
	end
	keVec =cellfun(@(Ke)Ke(tril(true(size(Ke)))),ke,'UniformOutput',0);