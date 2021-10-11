function [N,D] = ROM_Coefficient(Phi,nx,ny,dx,dy)
% - function for computing the Linear and non-linear coefficeint of the 
% - ROM equation
Ub   = Phi;
r = size(Phi,2);
[Ubx, Uby ] = parDiff(Ub,nx,ny,dx,dy);
D = (-Ubx'*Ubx.*dx- Uby'*Uby.*dy)./100;

u = Phi(1:nx*ny,1:r);
v = Phi(nx*ny+1:end,1:r);

ux = Ubx(1:nx*ny,1:r);
uy = Uby(1:nx*ny,1:r);

vx = Ubx(nx*ny+1:end,1:r);
vy = Uby(nx*ny+1:end,1:r);

N = cell(r,1);
for k = 1:r
    for i = 1:r
        for j = 1:r
            N{k}(j,i) = -sum(sum(dx.*u(:,k).*(u(:,i).*ux(:,j)+v(:,i).*uy(:,j)) +...
                dy.*v(:,k).*(u(:,i).*vx(:,j)+v(:,i).*vy(:,j))));
        end
    end
end



end
