function ret =  generateFoundationTri(refinements)

border = [0,-40,-40, 52,52,12,  12,    0; %Create boundary.
          0,  0,-30,-30, 0, 0,-0.5, -0.5]; %[X coordinates; Y
                                           %coordinates]


gd = [2,max(size(border)), border(1,:), border(2,:)]';
dl = decsg(gd);

%Create mesh
[p, e, t] = initmesh(dl);
for n = 1:refinements
  [p, e, t] = refinemesh(dl, p, e, t);
end

pdemesh(p,e,t);
xlim([min(border(1,:)) max(border(1,:))])
ylim([min(border(2,:)) max(border(2,:))])