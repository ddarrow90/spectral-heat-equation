function u = exactSolution(r,z,th,t)
w5 = [8.7714838159599540191228671334096,...
      12.338604197466943986082097644459,...
      15.700174079711671037587715595026,...
      18.980133875179921120770736748467,...
      22.217799896561267868824764947529,...
      25.430341154222704252022674430820];
rng(112374);      
q = randperm(length(w5));
c = 2*randperm(length(w5));
u = 0*r;

for j=1:length(w5)
    u = u + c(j)*exp(-((pi*q(j))^2 + w5(j)^2)*t)*sin(q(j)*pi*z).*besselj(5,w5(j)*r).*cos(5*th);
end

w7 = [11.086370019245083845762764435931,...
      14.821268727013171251365392327191,...
      18.287582832481726446143071367864,...
      21.641541019848400775128815411562,...
      24.934927887673022268865641506073];

q = randperm(length(w7));
c = 2*randperm(length(w7));

for j=1:length(w7)
    u = u + c(j)*exp(-((pi*q(j))^2 + w7(j)^2)*t)*sin(q(j)*pi*z).*besselj(7,w7(j)*r).*sin(7*th);
end

w16 = [21.085146113064718937846237094033,...
       25.417019006342758261561265988748];

q = randperm(length(w16));
c = 2*randperm(length(w16));

for j=1:length(w16)
    u = u + c(j)*exp(-((pi*q(j))^2 + w16(j)^2)*t)*sin(q(j)*pi*z).*besselj(16,w16(j)*r).*sin(16*th);
end

end
