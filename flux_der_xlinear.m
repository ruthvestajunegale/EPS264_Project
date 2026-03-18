function [ der ] = flux_der_xlinear(z,hSq,phi,psi)
global Te0 Drho step zf Dflux_der 
h = hSq^.5;
der = 1/2*h^3*psi*(((h^2-1)*(5+h^4*(phi-1)-phi+h^2*(5*phi-7)+...
    h^6*(3-5*phi+2*phi^2)))/(1+h^4*(phi-1))^2 -4*log(h))/(2*h);
if (z>0) && step
    der = der+Drho*Te0*(4*(h^2-1)*(h+h^3*(phi-1)))/(1+h^4*(phi-1))^2/(2*h);
end
if step < 0.5
    der = der+interp1(zf,Dflux_der,z);
end
end

