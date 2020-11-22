function [rho_app] = average_apparent_density(alpha_w, alpha_b)

% file:///C:/Users/benny/Downloads/Belinha,%20Jorge_%20Jorge,%20R.%20M.%20Natal_%20Reis%20Campos,%20J.%20C._%20Tavares,%20Jo%C3%A3o%20Manuel%20R.%20S._%20Vaz,%20M%C3%A1rio%20A.%20P%20-%20Biodental%20Engineering%20V_%20Proceedings%20of%20the%205th%20International%20Conference%20on%20Biodental%20Engineering%20(BIODENTAL%20201.pdf
%  2.1 g/cm3
rho_app_cortical = 2.1;

rho_app = (alpha_w / alpha_b) * rho_app_cortical;
end