%% Aqcuire Molulus distribution from Displacement Data by Adjoint Method
%%% Test the Plane Stress Problem of FGM under uniaxial stretch
clc; clear all
%%% - Add function lab -
addpath(genpath('./stenglib-master'));
addpath(genpath('./igfe-forward'));
addpath(genpath('./adjoint-method'));
%% Load Experiment Data
%%% Load ux and uy
load("U.mat","U");
%% Plane Stress Problem Descriptions
%%% Geometry & Mesh
geoL = 9; geoH = 9;
nelL = 100; nelH = 100;
[meshInfo,contourShow] = Mesh_Info(geoL,geoH,nelL,nelH);
%%% Dirichlet Condition
coord = meshInfo.coord;
fixEdge = find(coord(:,1)==0);
fixdof = sort([2*fixEdge-1;2*fixEdge(end)]); 
%%% Neunmann Condition
P = 0.1;
alldof = meshInfo.nDof; F = sparse(alldof,1);
loadPonit = find(coord(:,1)==geoL); loadDof = loadPonit*2-1;
FMag = P*meshInfo.eh ; F(loadDof,1) = FMag;
F(loadDof(1),1) = F(loadDof(1),1) - FMag/2;
F(loadDof(end),1) = F(loadDof(end),1) - FMag/2;
%%% Total Boundary Conditions
BCInfo.fixdof = fixdof;
BCInfo.F = F;
%% Initial Material Properties
%%% Initial Modulus Distribution Guess
iter_Evec = Modulus_Init_Gauss(meshInfo,6,1);
% iter_Evec = 1.4*ones(meshInfo.nNod,1);
%%% Poisson's Ratio
nu = 0.3;
materialInfo.nu = nu;
%%% Initial Update
materialInfo = Material_Update(materialInfo,iter_Evec);
% Assemble FEM Information
femInfo = FEM_Assemble(meshInfo,materialInfo,BCInfo);
% Solve the Forward Problem
forwU = Forward_Solver(femInfo);
% Get Lagarange Multipliers
lambda = Lambda_Distance(femInfo,forwU,U);
% Calculate the grad
grad = Get_Grad(femInfo,lambda,forwU);
%% BFGS Method
iter = 1;
max_iter = 1000;
E_min = 0.8;
E_max = 8;
H = eye(meshInfo.nNod,meshInfo.nNod);
figure;
Plot_Modulus(meshInfo,contourShow,iter_Evec);
while iter < max_iter
    % Calculate Gradient Direction
    dk = -H*grad;
    % Update Modulus Vector
    iter_Evec = iter_Evec+dk;
    minus_idx = find(iter_Evec<0);
    max_idx = find(iter_Evec>E_max);
    iter_Evec(minus_idx) = E_min;
    iter_Evec(max_idx) = E_max;
    iter = iter+1;
    materialInfo = Material_Update(materialInfo,iter_Evec,iter);
    Plot_Modulus(meshInfo,contourShow,iter_Evec);
    %%% Update Hessian Matrix
    % Assemble FEM Information
    femInfo = FEM_Assemble(meshInfo,materialInfo,BCInfo);
    % Solve the Forward Problem
    forwU = Forward_Solver(femInfo);
    % Get Lagarange Multipliers
    lambda = Lambda_Distance(femInfo,forwU,U);
    % Calculate the grad
    grad_next = Get_Grad(femInfo,lambda,forwU);
    if grad_next == 0
        break
    end
    y_k = grad_next-grad;
    % H = H+d_k*d_k'/(d_k'*y_k)-H*y_k*y_k'*H/(y_k'*H*y_k);
    H = H+(1+y_k'*H*y_k/((dk'*y_k)))*dk*dk'/(dk'*y_k)-(dk*y_k'*H+H*y_k*dk')/(dk'*y_k);
    grad = grad_next;
end