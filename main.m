clear; close all; clc;
rng(2); % 乱数シードの固定（再現性の確保）

% 波形を読み込み
[s1, fs] = audioread("s1.wav");
[s2, fs] = audioread("s2.wav");
s = [s1, s2].'; % sは音源×時間

% 波形をプロット
sigLen = size(s, 2); % 信号長
timeAx = 0:1/fs:(sigLen-1)/fs; % 時間軸
% figure; plot(timeAx, s(1, :)); grid on;
% figure; plot(timeAx, s(2, :)); grid on;

% 自分で混合
A = [0.8, -0.6;...
     0.7, 0.9]; % 混合行列
x = A*s; % 混合
% figure; plot(timeAx, x(1, :)); grid on;
% figure; plot(timeAx, x(2, :)); grid on;

% 再生
% sound(x(1, :), fs);

% ICA
iterNum = 30; % 反復回数
stepSize = 0.5; % 勾配降下法のステップサイズ（大きくすると不安定，小さくすると反復回数がたくさん必要）

% ICAの結果を使って分離
[W, J] = calc_ica(x, iterNum, stepSize);

for t = 1 : sigLen
    y(:, t) = W * x(:, t);
end

y = y/max(abs(y), [], 'all');

plot(J)

% % 波形表示
% figure; plot(timeAx, y(1, :)); grid on;
% figure; plot(timeAx, y(2, :)); grid on;
% 
% % 再生
% sound(y(1, :), fs);

% 保存
fileName = "./y1.wav";
audiowrite(fileName, y(1, :).', fs);
fileName = "./y2.wav";
audiowrite(fileName, y(2, :).', fs);