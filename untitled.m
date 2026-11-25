clc; clear; close all;
% Лабораторная работа № 2, вариант 13: f(x)=sinh(x), x in [-pi, pi].
% Ряд содержит постоянную составляющую и K пар гармоник.
N0 = 1024;
limit = 1;                           % допустимая СКО, %

%% Гармоника с целым числом периодов
kpInt = 2; Nint = 64;
[eRi, eCi] = fourierError(makeSignal('harm', Nint, kpInt), kpInt);
fprintf('Целое KP=%d, N=%d, K=%d: СКО действ. = %.12g%%, компл. = %.12g%%\n', ...
    kpInt, Nint, kpInt, eRi, eCi);

%% Гармоника с нецелым числом периодов (KP=2.5 по заданию)
kp = 2.5;
Kh = 3:30;
[ehR, ehC] = errorByK(makeSignal('harm', N0, kp), Kh);
KhOpt = Kh(find(ehR < limit, 1));
assert(~isempty(KhOpt), 'Увеличьте диапазон Kh.');
Nh = (4*KhOpt+1):80;                % строго N > 4K
ehNR = zeros(size(Nh)); ehNC = ehNR;
for q = 1:numel(Nh)
    [ehNR(q), ehNC(q)] = fourierError(makeSignal('harm', Nh(q), kp), KhOpt);
end
NhOpt = Nh(find(ehNR < limit, 1));

%% Функция варианта 13: sinh(x)
Ks = 1:floor((N0-1)/4);             % строго K < N0/4
[esR, esC] = errorByK(makeSignal('sinh', N0, 0), Ks);
KsOpt = Ks(find(esR < limit, 1));
assert(~isempty(KsOpt), 'При N=1024 порог 1%% не достигнут.');
Ns = (4*KsOpt+1):N0;
esNR = zeros(size(Ns)); esNC = esNR;
for q = 1:numel(Ns)
    [esNR(q), esNC(q)] = fourierError(makeSignal('sinh', Ns(q), 0), KsOpt);
end
NsOpt = Ns(find(esNR < limit, 1));

fprintf('\nНецелое KP=%.1f, N=%d: минимальное K=%d (СКО %.4f%%)\n', ...
    kp, N0, KhOpt, ehR(Kh == KhOpt));
fprintf('При K=%d: минимальное допустимое N с СКО < 1%%: %d (%.4f%%)\n', ...
    KhOpt, NhOpt, ehNR(Nh == NhOpt));
fprintf('\nsinh(x), N=%d: минимальное K=%d (СКО %.4f%%)\n', ...
    N0, KsOpt, esR(Ks == KsOpt));
fprintf('При K=%d: минимальное допустимое N с СКО < 1%%: %d (%.4f%%)\n', ...
    KsOpt, NsOpt, esNR(Ns == NsOpt));

disp('Гармоника KP=2.5: K, СКО действ., СКО компл. (%)');
disp([Kh(:), ehR(:), ehC(:)]);
disp('Гармоника KP=2.5: N, СКО действ., СКО компл. (%)');
disp([Nh(:), ehNR(:), ehNC(:)]);
disp('sinh(x): K, СКО действ., СКО компл. (%)');
disp([Ks(:), esR(:), esC(:)]);
disp('sinh(x): N, СКО действ., СКО компл. (%)');
disp([Ns(:), esNR(:), esNC(:)]);

drawError(Kh, ehR, ehC, limit, 'K', 'Гармоника KP=2.5, N=1024');
drawError(Nh, ehNR, ehNC, limit, 'N', sprintf('Гармоника KP=2.5, K=%d', KhOpt));
drawError(Ks, esR, esC, limit, 'K', 'sinh(x), N=1024');
drawError(Ns, esNR, esNC, limit, 'N', sprintf('sinh(x), K=%d', KsOpt));

[~, ~, y] = fourierError(makeSignal('sinh', N0, 0), KsOpt);
x = linspace(-pi, pi, N0+1);
figure('Name', 'Восстановление sinh(x)');
plot(x, sinh(x), 'k-', x, y, 'r--'); grid on;
xlabel('x'); ylabel('f(x)');
title(sprintf('sinh(x) и ряд Фурье, N=%d, K=%d', N0, KsOpt));
legend('Исходная функция', 'Восстановление');

function drawError(x, eR, eC, limit, xlabelText, titleText)
    figure('Name', titleText);
    plot(x, eR, 'b-', x, eC, 'r--', 'LineWidth', 1.4); hold on;
    yline(limit, 'k:'); grid on;
    xlabel(xlabelText); ylabel('СКО, %'); title(titleText);
    legend('Действительный ряд', 'Комплексный ряд', 'Порог 1 %');
end

function f = makeSignal(kind, N, kp)
    n = 0:N;                        % N уникальных отсчётов + правая граница
    if strcmp(kind, 'harm')
        f = sin(2*pi*kp*n/N);
    else
        f = sinh(-pi + 2*pi*n/N);
    end
end

function [errR, errC] = errorByK(f, Kvalues)
    errR = zeros(size(Kvalues)); errC = errR;
    for q = 1:numel(Kvalues)
        [errR(q), errC(q)] = fourierError(f, Kvalues(q));
    end
end

function [errR, errC, yR] = fourierError(f, K)
    N = numel(f)-1;
    assert(K < N/4, 'Требуется строгое условие K < N/4.');
    n = 0:N-1; m = 0:N;
    meanF = mean(f(1:N));           % крайняя точка f(N+1) не дублируется
    yR = meanF*ones(1,N+1);
    yC = yR;
    for k = 1:K
        angleN = 2*pi*k*n/N;
        a = 2*mean(f(1:N).*cos(angleN));
        b = 2*mean(f(1:N).*sin(angleN));
        c = mean(f(1:N).*exp(-1i*angleN));
        angleM = 2*pi*k*m/N;
        yR = yR + a*cos(angleM) + b*sin(angleM);
        yC = yC + 2*real(c*exp(1i*angleM));
    end
    scale = 100/(max(f)-min(f));
    errR = std((yR(2:N)-f(2:N))*scale); % оба края исключены
    errC = std((yC(2:N)-f(2:N))*scale);
end
function [errR, errC, yR] = fourierError(f, K)
    N = numel(f)-1;
    assert(K < N/4, 'Требуется строгое условие K < N/4.');
    n = 0:N-1; m = 0:N;
    meanF = mean(f(1:N));           % крайняя точка f(N+1) не дублируется
    yR = meanF*ones(1,N+1);
    yC = yR;
    for k = 1:K
        angleN = 2*pi*k*n/N;
        a = 2*mean(f(1:N).*cos(angleN));
        b = 2*mean(f(1:N).*sin(angleN));
        c = mean(f(1:N).*exp(-1i*angleN));
        angleM = 2*pi*k*m/N;
        yR = yR + a*cos(angleM) + b*sin(angleM);
        yC = yC + 2*real(c*exp(1i*angleM));
    end
    scale = 100/(max(f)-min(f));
    errR = std((yR(2:N)-f(2:N))*scale); % оба края исключены
    errC = std((yC(2:N)-f(2:N))*scale);
end
function [errR, errC, yR] = fourierError(f, K)
    N = numel(f)-1;
    assert(K < N/4, 'Требуется строгое условие K < N/4.');
    n = 0:N-1; m = 0:N;
    meanF = mean(f(1:N));           % крайняя точка f(N+1) не дублируется
    yR = meanF*ones(1,N+1);
    yC = yR;
    for k = 1:K
        angleN = 2*pi*k*n/N;
        a = 2*mean(f(1:N).*cos(angleN));
        b = 2*mean(f(1:N).*sin(angleN));
        c = mean(f(1:N).*exp(-1i*angleN));
        angleM = 2*pi*k*m/N;
        yR = yR + a*cos(angleM) + b*sin(angleM);
        yC = yC + 2*real(c*exp(1i*angleM));
    end
    scale = 100/(max(f)-min(f));
    errR = std((yR(2:N)-f(2:N))*scale); % оба края исключены
    errC = std((yC(2:N)-f(2:N))*scale);
end