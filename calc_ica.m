function [W, J] = calc_ica(x, iterNum, stepSize)
    sigLen = size(x, 2); % 信号長
    J = zeros(iterNum-1, 1);
    W = randn(2);
    I = eye(2); % Iをサイズ2×2の単位行列で定義

    for l = 1 : iterNum-1
        E = zeros(2); % 経験値の計算用
        for t = 1 : sigLen
            y(:, t) = W*x(:, t);
    %        p = y(:, t)./max(abs(y(:, t)), eps);
            p = tanh(y(:, t));
            R = p * y(:, t).';
            E = E + (1/sigLen)*R;
        end
        W = W - stepSize * (E - I) * W;
        J(l, 1) = -log(abs(det(W))) - (1/sigLen)*sum(log((1/pi)*sech(y)), 'all');
    end
end