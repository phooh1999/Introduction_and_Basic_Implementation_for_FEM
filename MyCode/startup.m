clc

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

addpath(genpath(strcat(pwd,sslash,'test_solver')));

fprintf('\n')
disp(' __  __       ______ ______ __  __  ')
disp(' |  \/  |     |  ____|  ____|  \/  |')
disp(' | \  / |_   _| |__  | |__  | \  / |')
disp(' | |\/| | | | |  __| |  __| | |\/| |')
disp(' | |  | | |_| | |    | |____| |  | |')
disp(' |_|  |_|\__, |_|    |______|_|  |_|')
disp('          __/ |                     ')
disp('         |___/                 MyFEM')
fprintf('\n\n')