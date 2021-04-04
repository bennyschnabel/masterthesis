clc; clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate fabric tensor based on mil tensor
% File name MIL tensor
fileNameMIL = 'Knochenprobe_2_06_600_M.dat';
%fileNameMIL = 'Knochenprobe_2_06_600_M_rot.csv';
% % File name for fabric tensor export file
fileNameFabric = 'Knochenprobe_2_06_600_H.csv';
%fileNameFabric = 'Knochenprobe_2_06_600_H_rot.csv';
% Calculate fabric tensors ['true', 'skip']
calcFabric = 'skip';

%% Calculate stiffness tensors based direct on fem input
% File name FEM data
fileNameFEM = 'analyze_Knochenprobe_2_0.6_0_eff_s.csv';
% Calc stiffness/ compliance tensors on direct FEM data input ['stiffness', 'compliance', 'skip']
calcCowin = 'skip';
% File name for Cowin coefficients export file
fileNameCowin = 'Knochenprobe_2_06_600_cowin_coefficients';
% File name for direct stiffness/ compliance tensors calculation based on Cowin coefficients
fileNameDirect = 'Knochenprobe_2_06_600_direct';

%% Calculate Cowin coefficients based on first x entrys of fem data
% File name for Cowin coefficients export
fileNameCowinCoeffTotal = 'Knochenprobe_2_06_600_cowin_coefficients_50_pinv.csv';
% Calc stiffness/ compliance tensors on direct FEM data input ['true','skip']
calcCowinFirstXFEM = 'true';
% Number of rows of FEM for calculation
numbRowsDM = 50;

% File name for direct stiffness/ compliance tensors calculation based on Cowin coefficients
fileNameEstimated = 'Knochenprobe_2_06_600_estimated';
% Calculate stiffness/ compliance tensor ['stiffness', 'compliance', 'skip']
calcTensor = 'skip';
% Mean or mode ['mean', 'mode', specificRow, total]
calcStatic = 'total';
% If specificRow specifie here
row = 1591;
% File name BVTVdata
fileNameBVTV = 'Knochenprobe_2_06_600_bvtv.csv';
% Threshold for BV/TV
%thresholdBVTV = 0.6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required packages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pkg load optim;
pkg load statistics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start to record calculation time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
c = cputime(); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fabric calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch (calcFabric)
  case {'true', 'TRUE'}
    printf ("Fabric calculation started\n");
    x = dlmread(fileNameMIL);
    %x = x(2:end,:);
    [sizeX] = size(x);
    
    fid = fopen(fileNameFabric, 'w+');
    string = 'DomainNo;H11;H12;H13;H21;H22;H23;H31;H32;H33;x11;x12;x13;x21;x22;x23;x31;x32;x33';
    fprintf(fid, string);
    fprintf(fid, '\n');
    fclose(fid);
    
    for i = 1 : 1 : sizeX(1)
      dispString = [num2str(i), '/', num2str(sizeX(1)) , ' - Fabric calculation'];
      disp(dispString)
      
      M = [[x(i,2), x(i,3), x(i,4)];...
      [x(i,5), x(i,6), x(i,7)]; ...
      [x(i,8), x(i,9), x(i,10)]];

      x11 = x(i,11);
      x12 = x(i,12);
      x13 = x(i,13);
      x21 = x(i,14);
      x22 = x(i,15);
      x23 = x(i,16);
      x31 = x(i,17);
      x32 = x(i,18);
      x33 = x(i,19);
      
      try (M);
        %disp('MIL tensor M is symmetric positive definite.')
        H = M^(-1/2);
        if (isnumeric(H) == 1)
          try chol(H);
            %disp('Fabric tensor H is symmetric positive definite.')
            string = [num2str(i - 1), ';', num2str(H(1,1)), ';', num2str(H(1,2)), ';', num2str(H(1,3)), ...
            ';', num2str(H(2,1)), ';', num2str(H(2,2)), ';', num2str(H(2,3)), ';', num2str(H(3,1)), ...
            ';', num2str(H(3,2)), ';', num2str(H(3,3)), ';', num2str(x11), ';', num2str(x12), ';', num2str(x13),...
            ';', num2str(x21), ';', num2str(x22), ';', num2str(x23), ';', num2str(x31), ';', num2str(x32), ...
            ';', num2str(x33)];
            fid = fopen(fileNameFabric, 'a+');
            fprintf(fid, string);
            fprintf(fid, '\n');
            fclose(fid);
          catch ME
            %disp('Fabric tensor H is not symmetric positive definite.')
            string = [num2str(i - 1), ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ...
            ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ...
            ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ...
            ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ...
            ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN)];
            fid = fopen(fileNameFabric, 'a+');
            fprintf(fid, string);
            fprintf(fid, '\n');
            fclose(fid);
          end
        else
          %disp('Fabric tensor H is complex')
          string = [num2str(i - 1), ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ...
          ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ...
          ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ...
          ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ...
          ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN)];
          fid = fopen(fileNameFabric, 'a+');
          fprintf(fid, string);
          fprintf(fid, '\n');
          fclose(fid);
        endif
      catch ME
        %disp('MIL tensor M is not symmetric positive definite')
        string = [num2str(i - 1), ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ...
        ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ...
        ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ...
        ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ...
        ';', num2str(NaN), ';', num2str(NaN), ';', num2str(NaN)];
        fid = fopen(fileNameFabric, 'a+');
        fprintf(fid, string);
        fprintf(fid, '\n');
        fclose(fid);
      end
    endfor
  otherwise
    printf ('Fabric calculation skipped\n');
endswitch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cowin coefficients calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch (calcCowin)
  case {'stiffness'}
    printf ('Cowin coefficients calculation started - stiffness\n');
    fileNameCowinS = strcat(fileNameCowin, '_stiffness.csv');
    fileNameDirect = strcat(fileNameDirect, '_stiffness.csv');
    
    x = dlmread(fileNameFabric);
    x = x(2:end,:);
    [sizeX] = size(x);
    
    y = dlmread(fileNameFEM, ',');
    y = y(2:end,:);
    y = sortrows(y, 1);
    
    fid = fopen(fileNameDirect, 'w+');
    string = 'DomainNo;C11;C21;C31;C41;C51;C61;C12;C22;C32;C42;C52;C62;C13;C23;C33;C43;C53;C63;C14;C24;C34;C44;C54;C64;C15;C25;C35;C45;C55;C65;C16;C26;C36;C46;C56;C66';
    fprintf(fid, string);
    fprintf(fid, '\n');
    fclose(fid);
    
    fid = fopen(fileNameCowinS, 'w+');
    string = 'DomainNo;k1;k2;k3;k4;k5;k6;k7;k8;k9';
    fprintf(fid, string);
    fprintf(fid, '\n');
    fclose(fid);

    for i = 1 : 1 : sizeX(1)
      string = [num2str(i), '/', num2str(sizeX(1)), ' - Cowin coefficients calculation'];
      disp(string)
      
      H = [[x(i,2), x(i,3), x(i,4)];...
      [x(i,5), x(i,6), x(i,7)]; ...
      [x(i,8), x(i,9), x(i,10)]];
      
      if ((H) == 0) == 1
        string = [num2str(i - 1), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN)];

        fid = fopen(fileNameDirect, 'a+');
        fprintf(fid, string);
        fprintf(fid, '\n');
        fclose(fid);
        
        string = [num2str(i - 1), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN)];
        
        fid = fopen(fileNameCowinS, 'a+');
        fprintf(fid, string);
        fprintf(fid, '\n');
        fclose(fid);
      elseif isnan(H(1,1)) == 1
        string = [num2str(i - 1), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN)];

        fid = fopen(fileNameDirect, 'a+');
        fprintf(fid, string);
        fprintf(fid, '\n');
        fclose(fid);
        
        string = [num2str(i - 1), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN)];
        
        fid = fopen(fileNameCowinS, 'a+');
        fprintf(fid, string);
        fprintf(fid, '\n');
        fclose(fid);
      else
        [V, lambda] = eig (H);
        lambda = [lambda(1,1); lambda(2,2); lambda(3,3)];
        lambda = sort(lambda, "descend");
        lambda1 = lambda(1) / (lambda(1) + lambda(2) + lambda(3));
        lambda2 = lambda(2) / (lambda(1) + lambda(2) + lambda(3));
        lambda3 = lambda(3) / (lambda(1) + lambda(2) + lambda(3));
        I2 = lambda1 * lambda2 + lambda1 * lambda3 + lambda2 * lambda3;
        
        if (lambda1 == lambda2)
          break;
        elseif (lambda2 == lambda3)
          break;
        elseif (lambda1 == lambda3)
          break;
        endif

        C = [[y(i,2), y(i,3), y(i,4), y(i,5), y(i,6), y(i,7)]; ...
        [y(i,8), y(i,9), y(i,10), y(i,11), y(i,12), y(i,13)]; ...
        [y(i,14), y(i,15), y(i,16), y(i,17), y(i,18), y(i,19)]; ...
        [y(i,20), y(i,21), y(i,22), y(i,23), y(i,24), y(i,25)]; ...
        [y(i,26), y(i,27), y(i,28), y(i,29), y(i,30), y(i,31)]; ...
        [y(i,32), y(i,33), y(i,34), y(i,35), y(i,36), y(i,37)]];
        
        CEntry = [C(1,1); C(2,2); C(3,3); C(1,2); C(1,3); C(2,3); C(6,6); C(5,5); C(4,4)];
        
        A = [[1, I2, 2 * lambda1, 2 * lambda1^2, lambda1^2, 2, 2 * I2, 4 * lambda1, 4 * lambda1^2]; ...
        [1, I2, 2 * lambda2, 2 * lambda2^2, lambda2^2, 2, 2 * I2, 4 * lambda2, 4 * lambda2^2]; ...
        [1, I2, 2 * lambda3, 2 * lambda3^2, lambda3^2, 2, 2 * I2, 4 * lambda3, 4 * lambda3^2]; ...
        [1, I2, lambda1 + lambda2, lambda1^2 + lambda2^2, lambda1 * lambda2, 0, 0, 0, 0]; ...
        [1, I2, lambda1 + lambda3, lambda1^2 + lambda3^2, lambda1 * lambda3, 0, 0, 0, 0]; ...
        [1, I2, lambda2 + lambda3, lambda2^2 + lambda3^2, lambda2 * lambda3, 0, 0, 0, 0]; ...
        [0, 0, 0, 0, 0, 1, I2, lambda1 + lambda2, lambda1^2 + lambda2^2]; ...
        [0, 0, 0, 0, 0, 1, I2, lambda1 + lambda3, lambda1^2 + lambda3^2]; ...
        [0, 0, 0, 0, 0, 1, I2, lambda2 + lambda3, lambda2^2 + lambda3^2]];
        
        cowinCoeff = pinv(A) * CEntry;
        %[cowinCoeff, sigma, r] = ols(CEntry, A);
        %[cowinCoeff] = pcr(A, CEntry, tol = 1e-6, maxit = 10000);
        %cowinCoeff = lsqnonneg(A, CEntry);
        %[cowinCoeff, flag, relres, iter, resvec] = bicgstab (A, CEntry);
        %cowinCoeff = gmres(A, CEntry, [], [], size(A,2));%;, [], [], [], [], z0);
        %[cowinCoeff, flag, relres, iter, resvec] = qmr(A, CEntry, [], 500, [], [], z0);
        
        CEntryCowin = A * cowinCoeff;
        
        CCowin = [[CEntryCowin(1), CEntryCowin(4), CEntryCowin(5), 0, 0, 0]; ...
        [CEntryCowin(4), CEntryCowin(2), CEntryCowin(6), 0, 0, 0]; ...
        [CEntryCowin(5), CEntryCowin(6), CEntryCowin(3), 0, 0, 0]; ...
        [0, 0, 0, CEntryCowin(9), 0, 0]; ...
        [0, 0, 0, 0, CEntryCowin(8), 0]; ...
        [0, 0, 0, 0, 0, CEntryCowin(7)]];
        
        string = [num2str(i - 1), ';', num2str(CCowin(1,1)), ';', num2str(CCowin(2,1)), ';', ...
        num2str(CCowin(3,1)), ';', num2str(CCowin(4,1)), ';', num2str(CCowin(5,1)), ';', ...
        num2str(CCowin(6,1)), ';', num2str(CCowin(1,2)), ';', num2str(CCowin(2,2)), ';', ...
        num2str(CCowin(3,2)), ';', num2str(CCowin(4,2)), ';', num2str(CCowin(5,2)), ';', ...
        num2str(CCowin(6,2)), ';', num2str(CCowin(1,3)), ';', num2str(CCowin(2,3)), ';', ...
        num2str(CCowin(3,3)), ';', num2str(CCowin(4,3)), ';', num2str(CCowin(5,3)), ';', ...
        num2str(CCowin(6,3)), ';', num2str(CCowin(1,4)), ';', num2str(CCowin(2,4)), ';', ...
        num2str(CCowin(3,4)), ';', num2str(CCowin(4,4)), ';', num2str(CCowin(5,4)), ';', ...
        num2str(CCowin(6,4)), ';', num2str(CCowin(1,5)), ';', num2str(CCowin(2,5)), ';', ...
        num2str(CCowin(3,5)), ';', num2str(CCowin(4,5)), ';', num2str(CCowin(5,5)), ';', ...
        num2str(CCowin(6,5)), ';', num2str(CCowin(1,6)), ';', num2str(CCowin(2,6)), ';', ...
        num2str(CCowin(3,6)), ';', num2str(CCowin(4,6)), ';', num2str(CCowin(5,6)), ';', ...
        num2str(CCowin(6,6))];
        
        fid = fopen(fileNameDirect, 'a+');
        fprintf(fid, string);
        fprintf(fid, '\n');
        fclose(fid);
        
        string = [num2str(i - 1), ';', num2str(cowinCoeff(1)), ';', num2str(cowinCoeff(2)), ';', ...
        num2str(cowinCoeff(3)), ';', num2str(cowinCoeff(4)), ';', num2str(cowinCoeff(5)), ';', ...
        num2str(cowinCoeff(6)), ';', num2str(cowinCoeff(7)), ';', num2str(cowinCoeff(8)), ';', ...
        num2str(cowinCoeff(9))];
        
        fid = fopen(fileNameCowinS, 'a+');
        fprintf(fid, string);
        fprintf(fid, '\n');
        fclose(fid);
      endif
    endfor
  case {'compliance'}
    printf ('Cowin coefficients calculation started - compliance\n');
    fileNameCowinS = strcat(fileNameCowin, '_compliance.csv');
    fileNameDirect = strcat(fileNameDirect, '_compliance.csv');
    
    x = dlmread(fileNameFabric);
    x = x(2:end,:);
    
    fid = fopen(fileNameDirect, 'w+');
    string = 'DomainNo;D11;D21;D31;D41;D51;D61;D12;D22;D32;D42;D52;D62;D13;D23;D33;D43;D53;D63;D14;D24;D34;D44;D54;D64;D15;D25;D35;D45;D55;D65;D16;D26;D36;D46;D56;D66';
    fprintf(fid, string);
    fprintf(fid, '\n');
    fclose(fid);
    
    id = fopen(fileNameCowinS, 'w+');
    string = 'DomainNo;k1;k2;k3;k4;k5;k6;k7;k8;k9';
    fprintf(fid, string);
    fprintf(fid, '\n');
    fclose(fid);
  otherwise
    printf ('Cowin coefficients calculation skipped\n');
endswitch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Cowin coefficients based on first x FEM entrys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch (calcCowinFirstXFEM)
  case {'true'}
    
    x = dlmread(fileNameFabric);
    x = x(2:end,:);
    [sizeX] = size(x);

    y = dlmread(fileNameFEM, ',');
    y = y(2:end,:);
    y = sortrows(y, 1);

    Ages = [];
    Cges = [];
    
    i = 1;
    j = 1;
    while j <= numbRowsDM
      string = [num2str(i), '/', num2str(sizeX(1)), ' - Cowin coefficients calculation'];
      disp(string)

      H = [[x(i,2), x(i,3), x(i,4)];...
      [x(i,5), x(i,6), x(i,7)]; ...
      [x(i,8), x(i,9), x(i,10)]];
      
      if ((H) == 0) == 1
        string = [num2str(i - 1), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN)];

      elseif isnan(H(1,1)) == 1
        string = [num2str(i - 1), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN)];
        
      else
        [V, lambda] = eig (H);
        lambda = [lambda(1,1); lambda(2,2); lambda(3,3)];
        lambda = sort(lambda, "descend");
        lambda1 = lambda(1) / (lambda(1) + lambda(2) + lambda(3));
        lambda2 = lambda(2) / (lambda(1) + lambda(2) + lambda(3));
        lambda3 = lambda(3) / (lambda(1) + lambda(2) + lambda(3));
        I2 = lambda1 * lambda2 + lambda1 * lambda3 + lambda2 * lambda3;
        
        if (lambda1 == lambda2)
          break;
        elseif (lambda2 == lambda3)
          break;
        elseif (lambda1 == lambda3)
          break;
        endif
        
        C = [[y(i,2), y(i,3), y(i,4), y(i,5), y(i,6), y(i,7)]; ...
        [y(i,8), y(i,9), y(i,10), y(i,11), y(i,12), y(i,13)]; ...
        [y(i,14), y(i,15), y(i,16), y(i,17), y(i,18), y(i,19)]; ...
        [y(i,20), y(i,21), y(i,22), y(i,23), y(i,24), y(i,25)]; ...
        [y(i,26), y(i,27), y(i,28), y(i,29), y(i,30), y(i,31)]; ...
        [y(i,32), y(i,33), y(i,34), y(i,35), y(i,36), y(i,37)]];

        CEntry = [C(1,1); C(2,2); C(3,3); C(1,2); C(1,3); C(2,3); C(6,6); C(5,5); C(4,4)];

        A = [[1, I2, 2 * lambda1, 2 * lambda1^2, lambda1^2, 2, 2 * I2, 4 * lambda1, 4 * lambda1^2]; ...
        [1, I2, 2 * lambda2, 2 * lambda2^2, lambda2^2, 2, 2 * I2, 4 * lambda2, 4 * lambda2^2]; ...
        [1, I2, 2 * lambda3, 2 * lambda3^2, lambda3^2, 2, 2 * I2, 4 * lambda3, 4 * lambda3^2]; ...
        [1, I2, lambda1 + lambda2, lambda1^2 + lambda2^2, lambda1 * lambda2, 0, 0, 0, 0]; ...
        [1, I2, lambda1 + lambda3, lambda1^2 + lambda3^2, lambda1 * lambda3, 0, 0, 0, 0]; ...
        [1, I2, lambda2 + lambda3, lambda2^2 + lambda3^2, lambda2 * lambda3, 0, 0, 0, 0]; ...
        [0, 0, 0, 0, 0, 1, I2, lambda1 + lambda2, lambda1^2 + lambda2^2]; ...
        [0, 0, 0, 0, 0, 1, I2, lambda1 + lambda3, lambda1^2 + lambda3^2]; ...
        [0, 0, 0, 0, 0, 1, I2, lambda2 + lambda3, lambda2^2 + lambda3^2]];
        
        Ages = [Ages; A];
        Cges = [Cges; CEntry];
        
        j = j + 1;

      endif

      i = i + 1;

    endwhile
  
    cowinCoeff = pinv(Ages) * Cges;

    fid = fopen(fileNameCowinCoeffTotal, 'w+');
    string = 'k1;k2;k3;k4;k5;k6;k7;k8;k9';
    fprintf(fid, string);
    fprintf(fid, '\n');
    fclose(fid);

    string = [num2str(cowinCoeff(1)), ';', num2str(cowinCoeff(2)), ';', ...
    num2str(cowinCoeff(3)), ';', num2str(cowinCoeff(4)), ';', num2str(cowinCoeff(5)), ';', ...
    num2str(cowinCoeff(6)), ';', num2str(cowinCoeff(7)), ';', num2str(cowinCoeff(8)), ';', ...
    num2str(cowinCoeff(9))];

    fid = fopen(fileNameCowinCoeffTotal, 'a+');
    fprintf(fid, string);
    fprintf(fid, '\n');
    fclose(fid);
  otherwise
    printf ('Cowin coefficients calculation skipped\n');
endswitch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate stiffness/ compliance tensor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch (calcTensor)
  case {'stiffness'}
    printf ('Stiffness tensor calculation started\n');

    fileNameEstimated = strcat(fileNameEstimated, '_stiffness.csv');
    fileNameCowinS = strcat(fileNameCowin, '_stiffness.csv');
    
    bvtv = dlmread(fileNameBVTV);
    bvtv = bvtv(2:end,:);
    
    cowinCoeff = dlmread(fileNameCowinS);
    cowinCoeff = cowinCoeff(2:end,:);
    
    [rowsBVTV] = find(bvtv(:,2) >= thresholdBVTV);
    rowsBVTV = [rowsBVTV];
    
    k1 = zeros(length(rowsBVTV), 1);
    k2 = zeros(length(rowsBVTV), 1);
    k3 = zeros(length(rowsBVTV), 1);
    k4 = zeros(length(rowsBVTV), 1);
    k5 = zeros(length(rowsBVTV), 1);
    k6 = zeros(length(rowsBVTV), 1);
    k7 = zeros(length(rowsBVTV), 1);
    k8 = zeros(length(rowsBVTV), 1);
    k9 = zeros(length(rowsBVTV), 1);
    
    for i = 1 : 1 : length(rowsBVTV)
      j = rowsBVTV(i);
      k1(i) = cowinCoeff(j,2);
      k2(i) = cowinCoeff(j,3);
      k3(i) = cowinCoeff(j,4);
      k4(i) = cowinCoeff(j,5);
      k5(i) = cowinCoeff(j,6);
      k6(i) = cowinCoeff(j,7);
      k7(i) = cowinCoeff(j,8);
      k8(i) = cowinCoeff(j,9);
      k9(i) = cowinCoeff(j,10);
    endfor
    
    k1mean = nanmean(k1);
    k1var = nanstd(k1, 0);
    k1mode = mode(k1);

    k2mean = nanmean(k2);
    k2var = nanstd(k2, 0);
    k2mode = mode(k2);

    k3mean = nanmean(k3);
    k3var = nanstd(k3, 0);
    k3mode = mode(k3);

    k4mean = nanmean(k4);
    k4var = nanstd(k4, 0);
    k4mode = mode(k4);

    k5mean = nanmean(k5);
    k5var = nanstd(k5, 0);
    k5mode = mode(k5);

    k6mean = nanmean(k6);
    k6var = nanstd(k6, 0);
    k6mode = mode(k6);

    k7mean = nanmean(k7);
    k7var = nanstd(k7, 0);
    k7mode = mode(k7);

    k8mean = nanmean(k8);
    k8var = nanstd(k8, 0);
    k8mode = mode(k8);

    k9mean = nanmean(k9);
    k9var = nanstd(k9, 0);
    k9mode = mode(k9);
        
    x = dlmread(fileNameFabric);
    x = x(2:end,:);
    [sizeX] = size(x);
    
    fid = fopen(fileNameEstimated, 'w+');
    string = 'DomainNo;C11;C21;C31;C41;C51;C61;C12;C22;C32;C42;C52;C62;C13;C23;C33;C43;C53;C63;C14;C24;C34;C44;C54;C64;C15;C25;C35;C45;C55;C65;C16;C26;C36;C46;C56;C66';
    fprintf(fid, string);
    fprintf(fid, '\n');
    fclose(fid);
    
    for i = 1 : 1 : sizeX(1)
      string = [num2str(i), '/', num2str(sizeX(1)), ' - Stiffness tensor calculation'];
      disp(string)
      
      H = [[x(i,2), x(i,3), x(i,4)];...
      [x(i,5), x(i,6), x(i,7)]; ...
      [x(i,8), x(i,9), x(i,10)]];
      
      if ((H) == 0) == 1
        string = [num2str(i - 1), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN)];

        fid = fopen(fileNameEstimated, 'a+');
        fprintf(fid, string);
        fprintf(fid, '\n');
        fclose(fid);
      elseif isnan(H(1,1)) == 1
        string = [num2str(i - 1), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN), ';', num2str(NaN), ';', num2str(NaN), ';', ...
        num2str(NaN)];

        fid = fopen(fileNameEstimated, 'a+');
        fprintf(fid, string);
        fprintf(fid, '\n');
        fclose(fid);
      else
        [V, lambda] = eig (H);
        lambda = [lambda(1,1); lambda(2,2); lambda(3,3)];
        lambda = sort(lambda, "descend");
        lambda1 = lambda(1) / (lambda(1) + lambda(2) + lambda(3));
        lambda2 = lambda(2) / (lambda(1) + lambda(2) + lambda(3));
        lambda3 = lambda(3) / (lambda(1) + lambda(2) + lambda(3));
        I2 = lambda1 * lambda2 + lambda1 * lambda3 + lambda2 * lambda3;
        
        if (lambda1 == lambda2)
          break;
        elseif (lambda2 == lambda3)
          break;
        elseif (lambda1 == lambda3)
          break;
        endif

        A = [[1, I2, 2 * lambda1, 2 * lambda1^2, lambda1^2, 2, 2 * I2, 4 * lambda1, 4 * lambda1^2]; ...
        [1, I2, 2 * lambda2, 2 * lambda2^2, lambda2^2, 2, 2 * I2, 4 * lambda2, 4 * lambda2^2]; ...
        [1, I2, 2 * lambda3, 2 * lambda3^2, lambda3^2, 2, 2 * I2, 4 * lambda3, 4 * lambda3^2]; ...
        [1, I2, lambda1 + lambda2, lambda1^2 + lambda2^2, lambda1 * lambda2, 0, 0, 0, 0]; ...
        [1, I2, lambda1 + lambda3, lambda1^2 + lambda3^2, lambda1 * lambda3, 0, 0, 0, 0]; ...
        [1, I2, lambda2 + lambda3, lambda2^2 + lambda3^2, lambda2 * lambda3, 0, 0, 0, 0]; ...
        [0, 0, 0, 0, 0, 1, I2, lambda1 + lambda2, lambda1^2 + lambda2^2]; ...
        [0, 0, 0, 0, 0, 1, I2, lambda1 + lambda3, lambda1^2 + lambda3^2]; ...
        [0, 0, 0, 0, 0, 1, I2, lambda2 + lambda3, lambda2^2 + lambda3^2]];
        
        switch (calcStatic)
          case {'mean'}
            cowinCoeffK = [k1mean; k2mean; k3mean; k4mean; k5mean; k6mean; k7mean; k8mean; k9mean];
          case {'mode'}   
            cowinCoeffK = [k1mode; k2mode; k3mode; k4mode; k5mode; k6mode; k7mode; k8mode; k9mode];
          case {'specificRow'}
            cowinCoeffK = [cowinCoeff(row,2); cowinCoeff(row,3); cowinCoeff(row,4); cowinCoeff(row,5); cowinCoeff(row,6); cowinCoeff(row,7); cowinCoeff(row,8); cowinCoeff(row,9); cowinCoeff(row,10)];
          case {'total'}
            cowinCoeffTotal = dlmread(fileNameCowinCoeffTotal);
            cowinCoeffTotal = cowinCoeffTotal(2:end,:);
            cowinCoeffK = [cowinCoeffTotal(1); cowinCoeffTotal(2); cowinCoeffTotal(3); cowinCoeffTotal(4); cowinCoeffTotal(5); cowinCoeffTotal(6); cowinCoeffTotal(7); cowinCoeffTotal(8); cowinCoeffTotal(9)];
          otherwise
            cowinCoeffK = [k1mean; k2mean; k3mean; k4mean; k5mean; k6mean; k7mean; k8mean; k9mean];
        endswitch 
        
        CEntry = A * cowinCoeffK;
        
        CEst = [[CEntry(1), CEntry(4), CEntry(5), 0, 0, 0]; ...
        [CEntry(4), CEntry(2), CEntry(6), 0, 0, 0]; ...
        [CEntry(5), CEntry(6), CEntry(3), 0, 0, 0]; ...
        [0, 0, 0, CEntry(9), 0, 0]; ...
        [0, 0, 0, 0, CEntry(8), 0]; ...
        [0, 0, 0, 0, 0, CEntry(7)]];
        
        string = [num2str(i - 1), ';', num2str(CEst(1,1)), ';', num2str(CEst(2,1)), ';', ...
        num2str(CEst(3,1)), ';', num2str(CEst(4,1)), ';', num2str(CEst(5,1)), ';', ...
        num2str(CEst(6,1)), ';', num2str(CEst(1,2)), ';', num2str(CEst(2,2)), ';', ...
        num2str(CEst(3,2)), ';', num2str(CEst(4,2)), ';', num2str(CEst(5,2)), ';', ...
        num2str(CEst(6,2)), ';', num2str(CEst(1,3)), ';', num2str(CEst(2,3)), ';', ...
        num2str(CEst(3,3)), ';', num2str(CEst(4,3)), ';', num2str(CEst(5,3)), ';', ...
        num2str(CEst(6,3)), ';', num2str(CEst(1,4)), ';', num2str(CEst(2,4)), ';', ...
        num2str(CEst(3,4)), ';', num2str(CEst(4,4)), ';', num2str(CEst(5,4)), ';', ...
        num2str(CEst(6,4)), ';', num2str(CEst(1,5)), ';', num2str(CEst(2,5)), ';', ...
        num2str(CEst(3,5)), ';', num2str(CEst(4,5)), ';', num2str(CEst(5,5)), ';', ...
        num2str(CEst(6,5)), ';', num2str(CEst(1,6)), ';', num2str(CEst(2,6)), ';', ...
        num2str(CEst(3,6)), ';', num2str(CEst(4,6)), ';', num2str(CEst(5,6)), ';', ...
        num2str(CEst(6,6))];
        
        fid = fopen(fileNameEstimated, 'a+');
        fprintf(fid, string);
        fprintf(fid, '\n');
        fclose(fid);
      endif
    endfor
  case {'compliance'}
    printf ('Compliance tensor calculation started\n');
    
    fileNameEstimated = strcat(fileNameEstimated, '_compliance.csv');
    fileNameCowin = strcat(fileNameCowin, '_compliance.csv');
  otherwise
    printf ('Stiffness/ compliance tensor calculation skipped\n');
endswitch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stop to record calculation time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = toc;
c = cputime() - c;

fprintf("Wall time was %f minutes \n", t / 60);
fprintf("CPU time was %f minutes \n", c / 60);