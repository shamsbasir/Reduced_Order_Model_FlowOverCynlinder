clear all, close all, clc
LW = 'linewidth'; FS = 'fontsize'; IN = 'interpret'; LT = 'latex';
FW = 'fontweight';B  = 'bold';
set(0,'defaulttextinterpreter','latex')
load CYLINDER.mat

%{
%% ------------- Animating the u velocity ------------
C1 = linspace(-.2,1.2,20);
C2 = linspace(-5,5,20);
figure(1)
for i=1:1:250
    plotCylinder(reshape(U(:,i),ny,nx),C1);
    colorbar
    title('U Velocity')
    hold off
    drawnow
end
%% ------------- Animating the y-gradeint of the u velocity ------------


figure(2)
for i=1:1:250
    [Ux,Uy] = gradient(reshape(U(:,i),ny,nx),dx,dy);
    [Vx,Vy] = gradient(reshape(V(:,i),ny,nx),dx,dy);
    VORT = Uy-Vx;

    plotCylinder(reshape(VORT,ny,nx),C2); hold on
    title('Vorticity')
    set(gca,'FontSize',15)
    drawnow
    hold off
end

%}
%%
A = [U;V];
h = dx*dy;
%% covariance matrix
C = A'*A*h;
%% POD modes
[Vc Sc VcT] = svd(C);
Sm = sqrt(Sc);
Phi = A*Vc*inv(Sm);

%% : V components of the first 6 dominant modes
Urow = size(U,1);
C1 = linspace(-.2,1.2,20);
C2 = linspace(-5,5,20);
figure
for i=1:1:6
    subplot(3,2,i)
    plotCylinder(reshape(Phi(Urow+1:end,i),ny,nx),C1);
    title(['POD mode [',num2str(i,'%1d'),']']);
    set(gca,FS,10)
    colorbar
end

%% : 25 eigenvalues of the covariance matrix
figure()
s = diag(Sc);
s = s(1:20);
semilogy(s,'-o',LW,1.6)
title("Singular values ");
xlabel("Singular value counts",FW,B,FS,12);
ylabel("Singular Values",FW,B,FS,12);

%% 2 :
dt = 0.2;
t = 20;
nsnap = t/(dt) +1;

for r =[3 5 7]
    y = Phi(:,1:r)'*A(:,nsnap)*h;
    u = Phi(1:Urow,1:r)*y;
    figure()
    subplot(1,2,1)
    % ROM
    plotCylinder(reshape(u(1:Urow,1),ny,nx),C1);
    colorbar
   title(['U using ROM with [',num2str(r,'%1d'),'] modes']);
    subplot(1,2,2)
colorbar
    % PDE
    plotCylinder(reshape(U(:,nsnap),ny,nx),C1);
    title(['Exact U at [',num2str(t,'%1d'),'s]']);
    colorbar
end

%% DMD
X1 = A(:,1:end-1);
X2 = A(:,2:end);
%% SVD and rank -r truncation
r = 21;
[U_dmd S_dmd V_dmd] = svd(X1,'econ');
Ur = U_dmd(:,1:r);
Sr = S_dmd(1:r,1:r);
Vr = V_dmd(:,1:r);

%% Build Atilde and DMD Modes
Atilde = Ur'*X2*Vr/Sr;
[W, eigs] = eig(Atilde);
Phi = X2*Vr/Sr*W; % DMD Modes

%% DMD Spectra
lambda = diag(eigs);                       % discrete-time DMD eigenvalues
omega = log(lambda)/dt;                    % the continuous-time DMD eigenvalues
[omega,I]=sort(omega,'descend','ComparisonMethod','real'); % sort the omegas
Phi= Phi(:,I(1:r));                        % sort the DMD Modes accordingly

%% 8 DMD Modes : 7 modes are requested but for the sake of subplots we plot 8
figure();
for i=1:1:8
    subplot(4,2,i);
    phi_re = real(Phi(:,i));
    plotCylinder(reshape(phi_re(Urow+1:end),ny,nx),C1);
    colorbar
    title(['DMD mode [',num2str(i,'%1d'),']']);
end


%% EigenValues
figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--',LW,1.2) % plot unit circle
hold on
h1 = scatter(real(lambda),imag(lambda),'ob','filled');
set(gca,FS,10)
axis equal


%% Compute DMD Solution
x1 = X1(:,1);
b = Phi\x1;
t = 0:dt:50;
time_dynamics = zeros(r,length(t));
for iter = 1:length(time_dynamics)
    time_dynamics (:,iter) = (b.*exp(omega*t(iter)));
end
X_dmd = Phi * time_dynamics;


%% DMD forcast
t = [45 49];
for i = 1:2
    figure()
    subplot(2,1,1)
    phi_r = real(X_dmd(:,t(i)/dt+1));
    plotCylinder(reshape(phi_r(1:Urow),ny,nx),C1);
    colorbar
    title(['POD forcast at Time = ',num2str(t(i)),' (s)'])
    subplot(2,1,2)
    plotCylinder(reshape(U(1:Urow,t(i)/dt+1),ny,nx),C1);
    title(['PDE at Time = ',num2str(t(i)),' (s)'])
    colorbar
end
