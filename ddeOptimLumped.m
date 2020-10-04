function Z1=ddeOptimLumped(X)
global nDays data sigma gamma1 tolN cf
beta=exp(X(1)); p=exp(-X(2)^2); tau=1+13*exp(-X(3)^2); in=exp(X(4));
%beta=X(1); pbar=X(2); tau=X(3); in=X(4);
nu=sigma+tau; pbar=p*exp(-gamma1*tau); 
lags=[sigma nu]; tspan=[cf nDays];
op=ddeset('reltol',10^(-tolN),'abstol',10^(-tolN));
sol=dde23(@ddeSet,lags,@hist,tspan,op);
Qin=sol.y(2,:); t=sol.x;
Z1=interp1(t,Qin,cf:nDays);

function z = ddeSet(~,y,Z)
ylag1=Z(1,1); ylag2=Z(1,2); 
z=zeros(2,1);
z(1)=pbar*exp(-beta*ylag2)-exp(-beta*ylag1)-gamma1*y(1)+1-pbar;
z(2)=beta*pbar*z(1)*exp(-beta*y(1));
end

    function z = hist(~)
        z=[in; data(cf)];
    end
end