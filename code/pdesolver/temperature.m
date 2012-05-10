function out = temperature()

PrepareTempInterpolation([6, 6;16,9]); %[Time, T; Time, T]
Tapr=[];
tStart=0;
dt=500;
tEnd=24*3600;
for t=tStart:dt:tEnd
	Tapr=[Tapr;Tout(t)];
end
PrepareTempInterpolation([6, -11;16,-5]);
Tdec=[];
for t=tStart:dt:tEnd
    Tdec=[Tdec;Tout(t)];
end
out=[Tapr , Tdec];

%Tapr = Tout(TidISekunder);
%Tdec = Tout(TidISekunder);

function ret = PrepareTempInterpolation(points)
global BreakTimes TemperatureSpline

%Points is a matrix that contains the known
%times in the first column and the assosciated temperatures
%in the second column. The function should be periodic on 24h.

splineCount = size(points,1)+1;
A = zeros(splineCount*2, splineCount*2);
b = zeros(splineCount*2,1);
BreakTimes = points(:,1);

%Prepare equation system

for n=1:(splineCount-1)
 A(2*n-1, [2*n-1, 2*n]) = [1, points(n,1)];
 A(2*n, [2*(n+1)-1, 2*(n+1)]) = [1, points(n,1)];
 b([2*n-1 2*n]) = points(n,2);
end

A(2*splineCount-1,[1,2,2*splineCount-1, 2*splineCount]) = ...
   [1, 0,-1,-24];

A(2*splineCount,[2, 2*splineCount]) = [1, -1];

%Solve for X.
TemperatureSpline = A\b;


function ret = Tout(Time)
global BreakTimes TemperatureSpline

t = mod(Time,24*3600)/3600;


foundSpline = 0;

for n = 1:(size(BreakTimes,1)-1)

 if((t > BreakTimes(n)) & (t < BreakTimes(n+1)))
   foundSpline = n+1;
 end
end

if(foundSpline == 0)

 if(t<= BreakTimes(1))
   foundSpline = 1;
 else
   foundSpline = size(BreakTimes,1)+1;
 end
end


ret = [1,t]*TemperatureSpline(2*foundSpline+[-1:0]);
