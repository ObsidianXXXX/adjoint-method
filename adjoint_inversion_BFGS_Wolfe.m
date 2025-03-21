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
% iter_Evec = 2*ones(meshInfo.nNod,1);
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
E_min = 0.5;
E_max = 8;
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
    [Evec_next,grad_next,lambda,forwU,iter] = Wolfe_Line_Search(U,forwU,grad,H,femInfo, ...
    materialInfo,E_min,E_max,alpha0,c1,c2);
    Plot_Modulus(meshInfo,contourShow,Evec_next);
    % Update Material Info
    materialInfo = Material_Update(materialInfo,Evec_next,iter);
    % Calculate cost function
    cost = {((U-forwU)'*(U-forwU))/2/norm_U};
    disp(cost);
    cost_his = [cost_his,cost];
    % Gradient Direction
    dk = Evec_next - iter_Evec;
    if grad_next == 0
        break
    end
    y_k = grad_next-grad;
    % Update Hessian Matrix
    % H = H+d_k*d_k'/(d_k'*y_k)-H*y_k*y_k'*H/(y_k'*H*y_k);
    H = H+(1+y_k'*H*y_k/((dk'*y_k)))*dk*dk'/(dk'*y_k)-(dk*y_k'*H+H*y_k*dk')/(dk'*y_k);
    grad = grad_next;
    iter_Evec = Evec_next;
end
figure;
plot(linspace(1,iter,iter),cell2mat(cost_his));