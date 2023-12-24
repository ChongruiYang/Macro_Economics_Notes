function y = gleg
    a = 0;
    b = 1;
    xgrid = [-0.906 -0.538 0 0.538 0.906];
    wgrid = [0.237 0.479 0.569 0.479 0.237];
    y = 0;
    for j = 1 : 5
        x = xgrid(j);
        z = (x+1)*(b-a)/2+a;
        y = y + (b-a)/2*wgrid(j)*f(z);
    end
end