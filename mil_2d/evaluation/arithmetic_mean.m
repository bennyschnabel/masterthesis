function A = arithmetic_mean(a)

xi = 0;
n = length(a);

for ii = 1 : n
    xi = xi + a(ii);
end

A = xi / n;

end