function [Fa] = flux_xlinear(z,hSq,phi,psi)
global Te0 Drho step zf Dflux 
h = hSq^.5;
Fa = -h^4*psi/8*(4*log(h)-...
    (h^2-1)*(4-phi+h^2*(3*phi-4))/(1+h^4*(phi-1)));
if (z>0) && step
    Fa = Fa+Drho*Te0*(h^2-1)^2/(1+h^4*(phi-1));
end
if step<0.5
    Fa = Fa+interp1(zf,Dflux,z);
end
end
