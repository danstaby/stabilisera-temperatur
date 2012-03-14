function ret = triderivative(pts, cor1, cor2, dim)
%Calculates the derivatives of linear triangular elements.
%cor1 and cor2 can be vectors


ret = 1/(pts(dim, cor1) - pts(dim, cor2));