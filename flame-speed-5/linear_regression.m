clear all
format long

function data = read_data ()
  data = dlmread("stat_points.csv", ";");
  [rows, cols] = size(data);
  params_num = 5;
  data = data(2:rows, 1:params_num);
endfunction

X = read_data();
X(:, 1) = log(X(:, 1));

[m, n] = size(X);
# n - число параметров, m - число векторов

b = X' * ones(m, 1);
A = X' * X;
c = A^-1 * b

Y = zeros(m, 1);
for i = 1 : m
  sum = 0;
  term = 0;
  for j = 1 : n
    term += c(j) * X(i, j);
  endfor
  sum += (term - 1) ^ 2;
  Y(i) = term;
endfor

Y