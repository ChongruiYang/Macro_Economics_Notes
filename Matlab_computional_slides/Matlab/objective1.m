function y = objective1(x, num)
    rng(1); % set random seed
    % num is number of simulations
    y = 0;
    for i = 1 : num
        z = rand;
        y = y + (z - x)^2;
    end
    y = y / num;
end