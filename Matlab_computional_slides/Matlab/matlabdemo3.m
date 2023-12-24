% Find the minimum, calculate integration, and find zero
p = fminsearch(@humps,.5)
Q = quadl(@humps,0,1)
z = fzero(@humps,.5)

% Calculate limits
syms x
limit((x^3 + 5)/(x^4 + 7))
limit(x^2 + 5, 3)

% Compute derivative
syms t
f = 3*t^2 + 2*t^(-2);
diff(f)

y = sin(x)
diff(y)

y = log10(x)
diff(y)

% Compute higher order derivative
f = x^2;
diff(f, 2)

% Finding the Maxima and Minima of a Curve
syms x
y = 2*x^3 + 3*x^2 - 12*x + 17;   % defining the function
ezplot(y, [-2, 2])
g = diff(y)
s = solve(g)
subs(y, 1), subs(y, -2)

% Compute integration
syms x 
int(2*x)
int(cos(x))
int(x, 4, 9)

f = x^2*cos(x);
ezplot(f, [-4,9])
a = int(f, -4, 9)
disp('Area: '), disp(double(a));