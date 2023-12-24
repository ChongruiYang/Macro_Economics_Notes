function y = midpoint(n)
    a = 0;
    b = 1;
    h = (b-a)/n;
    y = 0;
    for j = 1 : n
        x = a + (j - 1/2)*h;
        y = y + f(x)*h;
    end
end