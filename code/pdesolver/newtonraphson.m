function [uGuess wGuess TGuess res n]  = newtonraphson(fu, fw, fT, uGuess, ...
					   wGuess,TGuess, uFreeNodes, ...
					   wFreeNodes, TFreeNodes, ...
					   uRange,wRange,TRange,ITERMAX, epsilon)

res = epsilon+1;
n = 0;
resvec = zeros(max(size(wGuess)+size(uGuess)+size(TGuess)),1);
resnew = zeros(max(size(wGuess)+size(uGuess)+size(TGuess)),1);
x = zeros(max(size(resvec)),1);
x(uRange) = uGuess;
x(wRange) = wGuess;
x(TRange) = TGuess;
relaxation = 1;
MaxUpdate = 5;
wrongwrong = 0;
resvec(uRange,1) = LHSu(uGuess,wGuess,TGuess,uFreeNodes,wFreeNodes, ...
			TFreeNodes)- fu(uFreeNodes);
resvec(wRange,1) = LHSw(uGuess,wGuess,TGuess,uFreeNodes, wFreeNodes, ...
			  TFreeNodes) - fw(wFreeNodes);
resvec(TRange,1) = LHST(uGuess,wGuess,TGuess,uFreeNodes,wFreeNodes,TFreeNodes)...
        - fT(TFreeNodes);
res0 = 0.5*resvec'*resvec;
rvec = zeros(3,2);
dxnorm = [];
minAlpha = 1;

prevAlpha = 1;
prevG = 0;
prevRes = 0;
backtrackCount = 0;
iterCount = 0;

while ((sqrt(2*res) > epsilon || n < 2) && n < ITERMAX)
  res = 2*res0;
  MaxImprove = 13;
  relaxation = 1;
  m = 0;
  %if(wrongwrong == 1)
  %  disp('OK! Were here')
  %end
    
  J = calcJacobian(uGuess, wGuess, TGuess, uFreeNodes, wFreeNodes, ...
                     TFreeNodes, uRange, wRange, TRange);
 
  
  
  iterCount = iterCount + 1;
  
  dx = J\resvec;  
  xnew = x - dx;
    
  resnew(uRange,1) = LHSu(xnew(uRange),xnew(wRange),xnew(TRange),uFreeNodes, ...
			  wFreeNodes,TFreeNodes) - fu(uFreeNodes);
  resnew(wRange,1) = LHSw(xnew(uRange),xnew(wRange),xnew(TRange),uFreeNodes, wFreeNodes, ...
			TFreeNodes) - fw(wFreeNodes);
  resnew(TRange,1) = LHST(xnew(uRange),xnew(wRange),xnew(TRange),...
			  uFreeNodes,wFreeNodes,TFreeNodes)...
                          - fT(TFreeNodes);
  res = 0.5*resnew'*resnew;

  if(res < res0)
    %Do first jump
    x = xnew;
    resvec = resnew;
    uGuess = x(uRange);
    wGuess = x(wRange);
    TGuess = x(TRange);
    res0 = res;
    backtrackCount = 0;
    minAlpha = 1;
    prevAlpha = 1;
    prevG = res;
    iterCount = iterCount + 1;
    fprintf(1, '%i\t\t%6.4g\t%6.4g\t%6.4g\n', n, norm(resvec(TRange)),...
              norm(resvec(uRange)), norm(resvec(wRange)))

  else

    %Calculate derivative
    gprim = -resvec'*resvec;
    g1 = res;
    g0 = res0;

    if(backtrackCount == 0)
      %Minimize the quadratic function
      
      backtrackCount = backtrackCount+1;
      
      minAlpha = -gprim/(2*(g1-g0-gprim));
      
      %if(minAlpha > 0.5*prevAlpha)
	%minAlpha = 0.5*prevAlpha;
      %end

      
    else
      %Minimize the cubic function
      backtrackCount = backtrackCount+1;
      
      
      l1 = 1;
      l2 = prevAlpha;
      g1 = res;
      g2 = prevG;

      couf = 1/(l1-l2)*[1/l1^2, -1/l2^2;-l2/l1^2,l1/l2^2]*...
	     [g1-gprim*l1-res0;g2-gprim*l2-res0];

    

      minAlpha = gprim/(-couf(2) - sqrt(couf(2)^2-3*couf(1)*gprim));
      %alphavec = [0:0.01:1];
      %fundi = couf(1)*alphavec.^3+couf(2)*alphavec.^2+gprim* ...
	%      alphavec+g0;
      %figure(3)
      
      %plot(alphavec, fundi)
      %hold on
      %plot(0,g0,'*')
      %plot(l1, g1, '*')
      %plot(l2,g2, '*')
      %plot(minAlpha, couf(1)*minAlpha^3+couf(2)*minAlpha^2+gprim*minAlpha+g0,'o');
      %hold off
      %pause
      if(minAlpha < 0.1*prevAlpha(1,1))
	%minAlpha = 0.1*prevAlpha(1,1)
      elseif(minAlpha > 0.5*prevAlpha(1,1))
	minAlpha = 0.5*prevAlpha(1,1);
	
      end

      if(minAlpha < 2^(-15))
        minAlpha = 1;
        backtrackCount = 0;
      end

    end

    if(backtrackCount >= 2)
      
     if(mean(deltares)/jumpInitial < 1e-4)
         minAlpha = 1;
     end
    end

    %Perform jump
    iterCount = iterCount + 1;
    
    dx = J\resvec;
    xnew = x - minAlpha*dx;

    resnew(uRange,1) = LHSu(xnew(uRange),xnew(wRange),xnew(TRange), ...
			    uFreeNodes, ...            \
      wFreeNodes,TFreeNodes) - fu(uFreeNodes);
    resnew(wRange,1) = LHSw(xnew(uRange),xnew(wRange),xnew(TRange), ...
			    uFreeNodes, wFreeNodes, ...
                          TFreeNodes) - fw(wFreeNodes);
    resnew(TRange,1) = LHST(xnew(uRange),xnew(wRange),xnew(TRange),...
                          uFreeNodes,wFreeNodes,TFreeNodes)...
                          - fT(TFreeNodes);
    
    prevAlpha = minAlpha;
    x = xnew;
    res = 0.5*resnew'*resnew;
    prevG = res;
    resvec = resnew;
    uGuess = x(uRange);
    wGuess = x(wRange);
    TGuess = x(TRange);
    
    if(backtrackCount == 1)
      jumpInitial = -resnew'*resnew;
      deltares = [res0-res];
    else
      deltares = [deltares; res0-res];
    end
    res0 = res;
      fprintf(1, '%i\t\t%6.4g\t%6.4g\t%6.4g\tAlpha: %6.4g\n', n, norm(resvec(TRange)),...
              norm(resvec(uRange)), norm(resvec(wRange)),minAlpha)
    
  end

  
  n = n + 1;
end

