function var = variance(a, mu)

xi = 0;
n = length(a);
for ii = 1 : n
    xi = xi + (a(ii) - mu)^2;
end
var = xi / (n - 1);
end