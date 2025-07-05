root = fileparts (mfilename ('fullpath'));

here = pwd ();

opts.force_recompile = false;
opts.include_dir = fullfile (root, 'misc', 'include');

relpath = fullfile ('misc', 'pareto', 'private');
info = register_mex ([], relpath, 'stk_paretofind_mex', {}, {'pareto.h'});
info = register_mex (info, relpath, 'stk_isdominated_mex', {}, {'pareto.h'});
info = register_mex (info, relpath, 'stk_dominatedhv_mex', {'wfg.c'}, {'wfg.h'});
info = register_mex (info, relpath, 'stk_IWFG_mex', {'wfg.c','iwfg.c', 'pq.c'}, {'wfg.h', 'iwfg.h', 'pq.h'});

if exist ('OCTAVE_VERSION', 'builtin') == 5
    version_info = ['OCTAVE ' OCTAVE_VERSION];
else
    version_info = ['MATLAB ' version];
end

for k = 1:(length (info))
    gcc_version_warning = stk_init__compile ...
        (version_info, fullfile (root, info(k).relpath), ...
        opts, info(k).mexname, info(k).other_src, info(k).includes);
end




function gcc_version_warning = stk_init__compile ...
    (version_info, d, opts, mexname, other_src, includes)

mex_filename = [mexname '.' mexext];
mex_fullpath = fullfile (d, mex_filename);

src_filename = [mexname '.c'];

% List of all source files (not counting include files)
src_files = [{src_filename} other_src];

compile = opts.force_recompile;

if ~ compile
    
    % Compile if the MEX is not present or out-dated
    dir_mex = dir (mex_fullpath);
    compile = isempty (dir_mex);
    
    % Compile if different version of Matlab/Octave (better safe than sorry)
    try
        fid = fopen ([mex_fullpath '.info']);
        assert (strcmp (version_info, fgetl (fid)));
    catch
        compile = true;
    end
    
    if ~ compile
        
        % Compile if one of the source files is more recent than the MEX
        for k = 1:(length (src_files))
            % Look for src file in current directory
            dir_src = dir (fullfile (d, src_files{k}));
            if isempty (dir_src)
                error ('STK:stk_init__build_mex:FileNotFound', ...
                    sprintf ('Source file %s not found', src_files{k}));
            end
            if dir_mex.datenum < dir_src.datenum
                compile = true;
                break;
            end
        end
        
        if (~ compile) && (~ isempty (includes))
            
            % Compile if one of the include files is more recent than the MEX
            for k = 1:(length (includes))
                % Look for header file in current directory
                dir_hdr = dir (fullfile (d, includes{k}));
                if isempty (dir_hdr)
                    % Look for header file in include directory
                    dir_hdr = dir (fullfile (opts.include_dir, includes{k}));
                    if isempty (dir_hdr)
                        error ('STK:stk_init__build_mex:FileNotFound', ...
                            sprintf ('Header file %s not found', includes{k}));
                    end
                end
                if dir_mex.datenum < dir_hdr.datenum
                    compile = true;
                    break;
                end
            end
        end
    end
end

gcc_version_warning = [];
isoctave = (exist ('OCTAVE_VERSION', 'builtin') == 5);

persistent use_silent
if isoctave && (isempty (use_silent))
    use_silent = false;
end

if compile
    
    fprintf ('Compiling MEX-file %s... ', mexname);
    
    here = pwd ();
    
    try  % Safely change directory
        cd (d);
        
        include = sprintf ('-I%s', opts.include_dir);
        
        if isoctave
            mex (src_files{:}, include);
        else
            warning_state = warning ('off', 'MATLAB:mex:GccVersion_link');
            lastwarn ('');
            
            if isempty (use_silent)
                % Try with -silent (only supported in recent versions of Matlab)
                try
                    mex ('-silent', src_files{:}, include);
                    use_silent = true;
                catch
                    % Try without -silent
                    mex (src_files{:}, include);
                    use_silent = false;
                end
            elseif use_silent
                mex ('-silent', src_files{:}, include);
            else
                mex (src_files{:}, include);
            end
            
            [lastmsg, lastid] = lastwarn ();
            if ~ isempty (strfind (lastid, 'GccVersion'))  %#ok<STREMP> 
                gcc_version_warning.msg = lastmsg;
                gcc_version_warning.id = lastid;
            end
            warning (warning_state);
        end
        
        fid = fopen ([mex_fullpath '.info'], 'wt');
        fprintf (fid, version_info);
        fclose (fid);
        
        fprintf ('ok.\n');
        if isoctave
            fflush (stdout);
        end
        cd (here);
        
    catch
        cd (here);
        rethrow (lasterror ());
    end
end

end % function

function info = register_mex (info, relpath, mexname, other_src, includes)

if nargin < 4
    other_src = {};
end

if nargin < 5
    includes = {};
end

k = 1 + length (info);

info(k).relpath = relpath;
info(k).mexname = mexname;
info(k).other_src = other_src;
info(k).includes = [{'stk_mex.h'} includes];

end % function