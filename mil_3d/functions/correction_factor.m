% correction_factor
function [r] = correction_factor(numberOfOrientations)

%numberOfOrientations = [500; 1000];

excelFileName = 'result_bonej.xlsx';

% Evaluation custom implementation

mu_a = zeros(length(numberOfOrientations), 1);
mu_c = mu_a;
var_a = zeros(length(numberOfOrientations), 1);
var_c = var_a;

for jj = 1 : 1 : length(numberOfOrientations)

    folder = 'results';
    fileName = ['Knochenprobe2_1mm_1_', num2str(numberOfOrientations(jj)), '_*.*'];
    filePattern = fullfile(folder, fileName);
    files = dir(filePattern);

    a = zeros(length(files), 1);
    b = a;
    c = a;

    for ii = 1 : 1 : length(files)
        importData = table2array(readtable(files(ii).name));
        MIL = importData(:,1);
        theta = importData(:,2);
        phi = importData(:,3);

        x = zeros(length(MIL), 1);
        y = x;
        z = x;

        for kk = 1 : 1 : length(x)
            [x(kk), y(kk), z(kk)] = sc2cc(MIL(kk), theta(kk), phi(kk));
        end

        [~, radii, ~, ~, ~] = ellipsoid_fit( [ x y z ], '' );
        a(ii) = min(radii);
        c(ii) = max(radii);
    end

    mu_a(jj) = mittelwert(a);
    mu_c(jj) = mittelwert(c);
end

%% Get boneJ Information

[~,sheet_name]=xlsfinfo(excelFileName);

numberOfOrientationsBonej = [];
mu_a_bonej = [];
mu_c_bonej = [];
var_a_bonej = [];
var_c_bonej = [];

for k=1:numel(sheet_name)
    
    value = string(sheet_name{k});
    l = strfind(value,'_');
    b = extractBetween(value,1,l(1) - 1);
    numberOfOrientationsBonej = [numberOfOrientationsBonej, str2double(b)];
    %b = value(1:l(1));
    [~,~,data{k}] = xlsread('result_bonej.xlsx',sheet_name{k});
    
    dataSize = size(data{1, k});
    dataA = zeros(dataSize(1), 1);
    dataC = zeros(dataSize(1), 1);
    for ii = 1 : 1 : dataSize(1)
        a = str2double(data{1, k}{ii, 3});
        b = str2double(data{1, k}{ii, 4});
        c = str2double(data{1, k}{ii, 5});
        dataA(ii) = min([a; b; c]);
        dataC(ii) = max([a; b; c]);
    end
    
    mu_a_bonej = [mu_a_bonej, mittelwert(dataA)];
    mu_c_bonej = [mu_c_bonej, mittelwert(dataC)];
end

r_a = zeros(length(mu_a), 1);
r_c = r_a;

for ii = 1 : 1 : length(mu_a)
    r_a(ii) = mu_a_bonej(ii) / mu_a(ii);
    r_c(ii) = mu_c_bonej(ii) / mu_c(ii);
end

mu_r_a = mittelwert(r_a);
mu_r_c = mittelwert(r_c);

r = mittelwert([mu_r_a; mu_r_c]);

end