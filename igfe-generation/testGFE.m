%% Test the Plane Stress Problem of FGM under uniaxial stretch
clc; clear all
%%% - Add function lab -
addpath(genpath('../stenglib-master'));
%  		 - Geo & Mesh -
geoL = 9 ; geoH = 9; nelL = 40; nelH = 40; elementType = 'CPS4';
E1 =  2.0 ; E2 = 8.0; nu = 0.3 ; sigma = 1 ; N = 0.1 ;
Beta = (log(E2/E1)) / geoH ;
% Beta = (E2 - E1) / geoH ;
Eele = ones(nelL*nelH,1);
[meshInfo,contourShow] = Init_PDEsolver_Mesh_UniforRec_CPS4(geoL,geoH,nelL,nelH);
%  		 - Material Property-
materialInfo.type = 'LinearElastic2DPlaneStress';
materialInfo.E1 = E1 ;
materialInfo.E2 = E2 ;
materialInfo.nu = nu ;
materialInfo.Beta = Beta ;
%   	 - Boundary Condition-
coord = meshInfo.coord; eleSize= meshInfo.eleSize; numNod = meshInfo.nNod;
%   	 - Displacement -
fixEdge = find(coord(:,1)==0); fixdofs = sort([2*fixEdge-1;2*fixEdge(end)]);
%   	 - Loads -
alldof = meshInfo.nDof; F = sparse(alldof,1);
loadPonit = find(coord(:,1)==geoL); loadDof = loadPonit*2-1;
FMag = N*meshInfo.eh ; F(loadDof,1) = FMag;
F(loadDof(1),1) = F(loadDof(1),1) - FMag/2;
F(loadDof(end),1) = F(loadDof(end),1) - FMag/2;
%   	 - BC Assemble -
BCInfo.fixdof = fixdofs; BCInfo.F = F;
%  		 - Output Request -
OutputRequest.History{1} = 'Comp';
OutputRequest.FieldAtNode{1} = {{'At_Gauss'},{'Stress'}};
%  		 - Pre Information Assemble-
[GFEInfo] = GFE_Assemble(meshInfo,elementType,materialInfo,...
	OutputRequest);

state = PDESolver_FEA_MATLAB(GFEInfo,Eele,BCInfo);
U = state.U;
ux = reshape(state.U(1:2:end), nelL+1, nelH+1)';
uy = reshape(state.U(2:2:end), nelL+1, nelH+1)';
fx = reshape(state.RF(1:2:end), nelL+1, nelH+1)';
fy = reshape(state.RF(2:2:end), nelL+1, nelH+1)';
E_ext = E1+Beta*contourShow.PlotY;
E_ext = E1*exp(Beta*contourShow.PlotY);
cd ..\
save('ux.mat',"ux");
save('uy.mat',"uy");
save('fx.mat',"fx");
save('fy.mat',"fy");
save("U.mat","U");
save("E_ext.mat","E_ext");

