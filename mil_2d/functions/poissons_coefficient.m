function [nu] = poissons_coefficient(alpha_b, alpha_w)

% file:///C:/Users/benny/Downloads/marques2017.pdf

alpha_t = alpha_b + alpha_w;

nu = (0.0 * alpha_b + 0.3 * alpha_w) / alpha_t;

end