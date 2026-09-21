function root = setup_paths()
% Add the submission methods, figure generators, and bundled FINUFFT runtime.
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root,'aux_functions'), fullfile(root,'figures generation'));
if isfile(fullfile(root,'finufft',['finufft.' mexext]))
    addpath(fullfile(root,'finufft'));
end
if exist('finufft2d1','file') ~= 2 || exist('finufft2d2','file') ~= 2 || ...
        exist('finufft','file') ~= 3
    error(['A compatible FINUFFT MATLAB installation is required. ' ...
        'See README.md for setup instructions.']);
end
end
