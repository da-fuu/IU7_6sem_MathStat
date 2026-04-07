clear all;
x = [-2.54,-0.79,-4.27,-3.09,-3.82,-0.61,-0.64,-1.24,-1.73,-2.91,-1.48,-1.28,-0.37,-1.88,-2.19,-1.61,-1.52,-3.17,-1.36,-3.08,-3.11,-3.07,-1.57,-1.51,-2.37,-0.58,-3.05,-2.93,-1.01,-1.40,-2.06,-3.05,-1.84,-1.24,-1.89,-2.06,-1.59,-2.83,-1.07,-2.96,-3.17,-3.08,-0.49,-3.11,-3.14,-2.30,-3.99,-1.56,-1.28,-3.46,-2.63,-0.82,-2.18,-0.89,-3.08,-1.13,-1.62,-1.06,-2.98,-1.55,-1.49,-1.65,-1.45,-2.29,-0.85,-1.44,-2.87,-2.40,-2.13,-3.52,-1.42,-3.64,-3.47,-2.05,-2.39,-2.07,-0.80,-1.52,-3.92,-2.22,-0.78,-2.60,-1.78,-1.61,-1.65,-2.06,-3.33,-3.41,-1.97,-1.74,-2.04,0.01,-1.37,-3.15,-2.35,-3.66,-1.79,-2.56,-1.87,-1.06,-0.64,-2.49,-1.85,-1.40,-0.86,-0.17,-0.62,-2.85,-2.12,-1.17,-2.48,-1.65,-3.74,-2.87,-3.15,-1.89,-1.34,-4.33,-0.96,-1.79]; 
n = length(x);
Mmax = max(x);
Mmin = min(x);
R = Mmax - Mmin;
mu = sum(x) / n;
Ssq = sum((x - mu).^2) / (n - 1);
fprintf('Максимальное значение выборки: %.3f\n', Mmax);
fprintf('Минимальное значение выборки: %.3f\n', Mmin);
fprintf('Размах выборки: %.3f\n', R);
fprintf('Оценка математического ожидания: %.3f\n', mu);
fprintf('Оценка дисперсии: %.3f\n', Ssq);

m = floor(log2(n)) + 2;
[arr_n, arr_c] = hist(x, m);
fprintf('Количество промежутков: %d\n', m);
fprintf('|---|------------------|-----|\n');
fprintf('│ i │        J_i       │ n_i │\n');
fprintf('|---|------------------|-----|\n');
half_delta = R / (2 * m);
for i = 1 : length(arr_c)
    l = arr_c(i) - half_delta;
    r = arr_c(i) + half_delta;
    if (i ~= length(arr_c))
        fprintf('│%2d │ [%.3f, % .3f) │%4d │\n', i, l, r, arr_n(i));
    else
        fprintf('│%2d │ [%.3f, % .3f] │%4d │\n', i, l, r, arr_n(i));
    end    
end
fprintf('|---|------------------|-----|\n');

sigma = sqrt(Ssq);
points = linspace(mu - 5*sigma, mu + 5*sigma, 500);
f = normpdf(points, mu, sigma);
norm_arr_n = arr_n / (sum(arr_n) * half_delta * 2);
figure('Name', 'Гистограмма и функция плотности')
bar(arr_c, norm_arr_n, 1);
hold on;
plot(points, f, 'LineWidth', 2, 'color', 'red');
xlabel('Значение выборки');
ylabel('Плотность');
legend('Гистограмма', 'N(\mu, \sigma^2)', 'Location', 'northwest');
hold off;

F = normcdf(points, mu, sigma);
[Fn, xn] = ecdf(x);
figure('Name', 'Функции распределения')
stairs(xn, Fn, 'LineWidth', 2); 
hold on;
plot(points, F, 'LineWidth', 1, 'color', 'red'); 
xlabel('Значения выборки');
ylabel('Вероятность');
legend('Эмпирическая', 'N(\mu, \sigma^2)', 'Location', 'northwest');
hold off;
