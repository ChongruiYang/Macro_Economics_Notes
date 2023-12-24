function y = objective2(w, num)
    rng(1); % set random seed
    % num is number of simulations
    mu = 1.06;
    sigma = 0.2;
    R = 1.01;
    y = 0;
    for i = 1 : num
        Z = normrnd(mu,sigma);
        y = y + exp(-((1-w)*R + w*Z));
    end
    y = y / num;
end