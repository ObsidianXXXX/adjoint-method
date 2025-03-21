function eleCen = PDEInfo_ELeCen(meshInfo)
	%%%===========================Copyright==================================%%%
	%%%   First Version Feb. 2020
	%%%
	%%%   Xueyan Hu <huxueyan@zju.edu.cn>
	%%%   Institute of Applied Mechanics,Zhejiang University
	%%%
	%%%===========================Description================================%%%
	%%% This function is to calculate the center coordinates for each element
	%%% ----INPUT
	%%%  nods_coor: <vector,[Num of Nodes , 2]> = coordinates of each node
	%%%  EleNodesID: <vector>,[Num of element , 4]> = node ID for each element
	%%% ----OUTPUT
	%%%  eleCen: <vector>,[Num of element , 2]> = center coordinate for each
	%%%  									     element
	%%%
	%%%======================================================================%%%
	nods_coor = meshInfo.coord; EleNodesID = meshInfo.eleNodsID;
	eleCen = zeros(length(EleNodesID(:,1)),2);
	for m = 1:length(EleNodesID(:,1))
		nodeId = EleNodesID(m,:);
		eleCen(m,1) = mean(nods_coor(nodeId,1));
		eleCen(m,2) = mean(nods_coor(nodeId,2));
	end