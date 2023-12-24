
% Use a for-loop to calculate 1+3+5+7+9
sum1 = 0;
for k = 1:2:9
 sum1 = sum1+k;
end
sum1

% Determine whether a given year is a leap year 
% try to change the given value of nyear and observe the outcome
nyear = 1975;
if (mod(nyear, 400) == 0)
 fprintf('%6u is a leap year', nyear)
elseif (mod(nyear,4) == 0) && (mod(nyear,100) ~= 0) 
 fprintf('%6u is a leap year', nyear)
else
 fprintf('%6u is not a leap year', nyear)
end

% Find N such that 3^N > 100 and 3^(N-1) < 100 
x = 3;
iter = 1;
while (x < 100)
 x = x*3
 iter = iter + 1
end
iter

% Construct a function to calculate 2*x^2 + 3*x + 7 for x = 1, 2, ..., 5
for n = 1:5
 x = n*0.1;
 z = myfunc1(x);
 fprintf('x = %4.2f f(x) = %8.4f \r',x,z)
end
