function ret = reflections(TTL)
global WindowCount Intensity a

a = 0.1;
k = 1/(1-a^2);
WindowCount = 3;

Intensity = 0;
doWindow(1-a, 1, 0, TTL);

disp(['Numerical result: ' num2str(Intensity)])
disp(['Analytical result: ' num2str((1-a)^3*k^2/(1-k^2*a^2*(1-a)^2))])

function doWindow(Iin, wID, lwID, TTL)
global WindowCount Intensity a

if(TTL == 0)
  return
end
if(wID+1 == lwID)
  k = a;
else
  k = 1-a;
end


if(wID ~= WindowCount)
  doWindow(k*Iin, wID+1, wID, TTL-1);
else
  Intensity = Intensity + Iin;
end

if(wID > 1)
  doWindow((1-k)*Iin, wID-1, wID, TTL-1); 
end
