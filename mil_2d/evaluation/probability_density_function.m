function [phi] = probability_density_function(x, sigma, mu)

phi = 1 / (sigma * sqrt(2 * pi)) * exp(-0.5 * ((x - mu) / sigma).^2);

end