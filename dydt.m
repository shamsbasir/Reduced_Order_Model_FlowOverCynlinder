function dy = dydt(y,Q,D,r)
% -- ODE 

dy = zeros(r,1);

for i = 1:r
    dy(i) = y'*Q{i}*y +D(i,:)*y;
end

end
