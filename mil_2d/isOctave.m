function retval = isOctave
% https://octave.org/doc/v4.0.1/How-to-distinguish-between-Octave-and-Matlab_003f.html
  persistent cacheval;  % speeds up repeated calls

  if isempty (cacheval)
    cacheval = (exist ("OCTAVE_VERSION", "builtin") > 0);
  end

  retval = cacheval;
end