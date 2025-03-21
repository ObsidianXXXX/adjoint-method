function grad = Get_Grad(femInfo,lambda,forwU)
%%%===========================Copyright==================================%%%
	%%%   Version Mar. 2024
	%%%
	%%%   Gengxuan Zhu <12324117@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This function is to calculate the grad of objective function
	%%%======================================================================%%%
	%%%======================================================================%%%
    meshInfo = femInfo.meshInfo;
    eleDofID = meshInfo.eleDofID;
    IDlist = reshape(meshInfo.eleNodsID',[],1);
    lambdaEle = num2cell(lambda(eleDofID),2);
    forwUEle = num2cell(forwU(eleDofID),2);
    Bgauss = femInfo.Bgauss;
    detJgauss = femInfo.detJgauss;
    D0 = femInfo.D0;
    [N,~,~,w] = shapeFunAtGauss();
    %% Calculate dKdE
    for i = 1:length(w)
		Bgauss_i = Bgauss(:,i);
        detJgauss_i = detJgauss(:,i);
        BD0BJgauss_i = cellfun(@(B,detJ)B'*D0*B*detJ,Bgauss_i,...
			detJgauss_i,'UniformOutput',0);
	    wi = w(i);
        N_i = N(i,:);
		if i==1 
            dKdE1 = cellfun(@(BD0BJ)wi*BD0BJ*N_i(1),BD0BJgauss_i,'UniformOutput',0);
            dKdE2 = cellfun(@(BD0BJ)wi*BD0BJ*N_i(2),BD0BJgauss_i,'UniformOutput',0);
            dKdE3 = cellfun(@(BD0BJ)wi*BD0BJ*N_i(3),BD0BJgauss_i,'UniformOutput',0);
            dKdE4 = cellfun(@(BD0BJ)wi*BD0BJ*N_i(4),BD0BJgauss_i,'UniformOutput',0);    
        else
			dKdE1 = cellfun(@(dKdE1,BDBJ)dKdE1+wi*BDBJ*N_i(1),dKdE1,BD0BJgauss_i,'UniformOutput',0);
            dKdE2 = cellfun(@(dKdE2,BDBJ)dKdE2+wi*BDBJ*N_i(2),dKdE2,BD0BJgauss_i,'UniformOutput',0);
            dKdE3 = cellfun(@(dKdE3,BDBJ)dKdE3+wi*BDBJ*N_i(3),dKdE3,BD0BJgauss_i,'UniformOutput',0);
            dKdE4 = cellfun(@(dKdE4,BDBJ)dKdE4+wi*BDBJ*N_i(4),dKdE4,BD0BJgauss_i,'UniformOutput',0);
        end   
    end
    gradEle = ...
    [cellfun(@(dKdE1,lambdaEle,forwUEle)lambdaEle*dKdE1*forwUEle',dKdE1,lambdaEle,forwUEle,'UniformOutput',1),...
     cellfun(@(dKdE2,lambdaEle,forwUEle)lambdaEle*dKdE2*forwUEle',dKdE2,lambdaEle,forwUEle,'UniformOutput',1),...
     cellfun(@(dKdE3,lambdaEle,forwUEle)lambdaEle*dKdE3*forwUEle',dKdE3,lambdaEle,forwUEle,'UniformOutput',1),...
     cellfun(@(dKdE4,lambdaEle,forwUEle)lambdaEle*dKdE4*forwUEle',dKdE4,lambdaEle,forwUEle,'UniformOutput',1)];
    
    grad = accumarray(IDlist,reshape(gradEle',[],1));
end