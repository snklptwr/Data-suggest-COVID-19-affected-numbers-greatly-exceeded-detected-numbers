function Z1=ddeOptimContinuum(X)
global tolN data nDays sigma tau gamma1 cf 
a=exp(X(1)); m=2.1+8*exp(-X(2)^2); p=exp(-X(3)^2); tau=1+13*exp(-X(4)^2); in=exp(X(5)); 
%a=X(1); m=X(2); pbar=X(3); tau=X(4); in=X(5); 
nu=sigma+tau; pbar=p*exp(-gamma1*tau);
lags=[sigma nu]; tspan=[cf nDays];
op=ddeset('reltol',10^(-tolN),'abstol',10^(-tolN));
sol=dde23(@ddeSet,lags,@hist,tspan,op);
P=sol.y(2,:); t=sol.x;
Z1=interp1(t,P,cf:nDays);

function z = ddeSet(~,y,Z)
ylag1=Z(1,1); ylag2=Z(1,2);
z=zeros(2,1); 
    function zG=Gfun(yArg)
             %zG1=-cg * (exp(cg * yArg) * cg ^ 2 * yArg ^ 2 * expint(cg * yArg) + 0.2e1 * exp(cg * yArg) * cg * yArg * expint(cg * yArg) - cg * yArg - 0.1e1);
             %if yArg<0; disp('oh yes'); end
             g2A=gamma_incomplete(abs(a*yArg),2-m);
             zG=yArg^(m-1)*a^m*exp(a*yArg)*g2A+yArg^(-2+m)*a^(m-1)*exp(a*yArg)*g2A*m-a^(m-1)*yArg^(-2+m)*exp(a*yArg)*g2A-a;
    end

z(1)=-Gfun(ylag1)+pbar*Gfun(ylag2)-gamma1*y(1)+(1-pbar)*a/(m-2);
z(2)=pbar*z(1)*Gfun(y(1));
end

  function z = hist(~)
      z=[10^(-2)*in; data(cf)];
  end
end