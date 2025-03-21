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
nelL = 2; nelH = 3;
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
%%% Poisson's Ratio
nu = 0.3;
materialInfo.nu = nu;
E_min = 1e-02;
iter_rho = Rho_Init_Gauss(meshInfo,2,-2);
rhoInfo.min = E_min;
%%% Initial Update
rhoInfo = Rho_Update(rhoInfo,iter_rho);
iter_Evec = E_min+iter_rho.^2;
materialInfo = Material_Update(materialInfo,iter_Evec);
% Assemble FEM Information
femInfo = FEM_Assemble(meshInfo,materialInfo,BCInfo);
% Solve the Forward Problem
forwU = Forward_Solver(femInfo);
% Get Lagarange Multipliers
lambda = Lambda_Distance(femInfo,forwU,U);
% Calculate the grad and conjugate direction
grad = 2*Get_Grad(femInfo,lambda,forwU).*iter_rho;
dk = -grad;
%% BFGS Method
iter = 1;
max_iter = 10000;
eps = 1e-07;
c1 = 0.0001;
c2 = 0.9;
alpha0 = 1;
H = eye(meshInfo.nNod,meshInfo.nNod);
norm_U = U'*U;
cost = {((U-forwU)'*(U-forwU))/2/norm_U};
cost_his = [cost];
disp(cost);
figure;
Plot_Modulus(meshInfo,contourShow,iter_Evec);
while iter < max_iter
    [rho_next,Evec_next,grad_next,forwU,iter] = CG_Line_Search(U,forwU,grad,dk,femInfo, ...
    materialInfo,rhoInfo,alpha0,c1,c2);
    Plot_Modulus(meshInfo,contourShow,Evec_next);
    % Update Material Info
    rhoInfo = Rho_Update(rhoInfo,rho_next,iter);
    materialInfo = Material_Update(materialInfo,Evec_next,iter);
    % Calculate cost function
    cost = {((U-forwU)'*(U-forwU))/2/norm_U};
    disp(cost);
    cost_his = [cost_his,cost];
    if cell2mat(cost) < eps
        break
    end
    % Update Conjugate Direction
    beta = grad_next'*grad_next/(grad'*grad);
    dk = -grad_next+beta*grad;
    grad = grad_next;
end
figure;
plot(linspace(1,iter,iter),cell2mat(cost_his));