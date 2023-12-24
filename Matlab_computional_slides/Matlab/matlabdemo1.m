% cd 'C:\Users\b134468\Dropbox\Course\ECON5170_Computational_methods\Matlab\Lecture 1'
cd 'C:\Users\EconUser\Dropbox\Course\ECON5170_Computational_methods\Matlab\Lecture 1'

% Load magic_matrix data
load magic_matrix.txt
A = magic_matrix;
% Calculate the summation of the magic matrix by row and by column
sum(A)
sum(A')';
% Replace the element in the second row and third column of the magic matrix by zero
A(2,3) = 0;
% Remove the second row of the magic matrix
A(2,:) = [];

% Generate a 2X2 matrix of random normal variables that follow N(1, 2^2)
rng(20);
B = (randn(2,2)*2+1);

% Calculate the determinant, eigenvalues, and inverse of the matrix you just generated
det(B)
inv(B)
eig(B)
% 7.	Calculate the mean and standard deviation of each row
mu = mean(B');
sigma = std(B');

% Generate a vector n that contains elements from 1 to 5. Build a table of n, n^2, and ln(n)
n = (1:5)';
C = [n n.^2 log(n)];
% Find prime numbers in the magic matrix
k = find(isprime(A))';
A(k)

% Draw sin(x), sin(x-0.25), and sin(x-0.5), X in [0, 2pi] on the same graph, use different color strings and marker styles to distinguish them, and add legend
x = 0:pi/10:2*pi;
y = sin(x);
y2 = sin(x-.25);
y3 = sin(x-.5);
plot(x,y,'r:+',x,y2,'b:o',x,y3,'g:x')
legend('sin(x)','sin(x-.25)','sin(x-.5)')

% Draw mesh and surface plots for Z = 1/sqrt(X^2 + Y^2), X in [-8, 8] and Y in [-8, 8]
[X,Y] = meshgrid(-8:.5:8); 
Z = 1./sqrt(X.^2 + Y.^2);
mesh(X,Y,Z,'EdgeColor','black')

surf(X,Y,Z)
colormap hsv
colorbar


[x,y,z] = peaks;
pcolor(x,y,z)
shading interp
hold on
contour(x,y,z,20,'k')
hold off