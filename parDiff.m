function [Ubx, Uby ]= parDiff(Ub,nx,ny,dx,dy)
% -- function for taking partial x and partial y derivatives

nt  = size(Ub,2);
np = nx*ny;
% allocate arrays
Ubx = 0.*Ub;
Uby = 0.*Ub;

for i = 1:nt
    % separate the u and v components
    Ub_u = Ub(1:np,i) ;
    Ub_v = Ub(np+1:end,i);
    % take the derivative
    [Ub_ux,Ub_uy]= gradient(reshape(Ub_u,ny,nx),dx,dy);
    
    Ubux = reshape(Ub_ux,ny*nx,1);
    Ubuy = reshape(Ub_uy,ny*nx,1);
    
    [Ub_vx,Ub_vy]= gradient(reshape(Ub_v,ny,nx),dx,dy);
    Ubvx = reshape(Ub_vx,ny*nx,1);
    Ubvy = reshape(Ub_vy,ny*nx,1);
    
    Ubx(1:np,i)     = Ubux;
    Ubx(np+1:end,i) = Ubvx;
    
    Uby(1:np,i)     = Ubuy ;
    Uby(np+1:end,i) = Ubvy;
end


end
