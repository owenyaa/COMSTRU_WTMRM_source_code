function x=l1_solution(SSS,dR,mu)
% initial solution
[U,ss,V]=svd(SSS);
x=V(:,1)*(U(:,1)'*dR)/ss(1,1);
eps=1e-8;
for itr=1:1000
    xtemp=x;
    W=mu*diag(1./(abs(x)+eps*max(abs(x))));
    x=(SSS'*SSS+W)\(SSS'*dR);
    if norm(x-xtemp)/norm(x)<1e-10
        break;
    end
end