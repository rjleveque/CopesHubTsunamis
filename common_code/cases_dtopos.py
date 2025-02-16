

def make_all_cases_dtopos(dtopo_dir, dtopo_files, xgeoclaw_path, runs_dir='.'):
    """
    Output: *caselist*, a list of cases to be run.
    Each case should be dictionary of any parameters needed to set up an
    individual case.  These will be used by run_one_case_dtopo.

    For this example, each dtopo_file in dtopo_files corresponds to a dtopo
    file.  GeoClaw will be run for each of these earthquake sources.
    A unique directory will be created for each run, with names
    based on dtopo_file, and residing within runs_dir.
    The _output and _plots directories will be within the run directory.

    """
    import os

    # Create a list of the cases to be run:
    caselist = []

    for dtopofile in dtopo_files:
        dtopo_name = os.path.splitext(os.path.split(dtopofile)[-1])[0]

        # Create all directories needed first, in case there's a problem:

        runs_dir = os.path.abspath(runs_dir)
        outdir = os.path.join(runs_dir, 'geoclaw_outputs/_output_%s' \
                    % dtopo_name)
        plotdir = os.path.join(runs_dir, 'geoclaw_plots/_plots_%s' \
                    % dtopo_name)
        os.system('mkdir -p %s' % outdir)
        print('Created %s' % outdir)
        os.system('mkdir -p %s' % plotdir)
        print('Created %s' % plotdir)

        case = {}

        case['case_name'] = dtopo_name
        case['outdir'] = outdir

        #case['xclawcmd'] = None  # if None, will not run code
        case['xclawcmd'] = xgeoclaw_path  # executable created by 'make .exe'

        # setrun parameters:
        case['setrun_file'] = 'setrun_case.py'
        case['dtopofiles'] = [[3, dtopofile]]

        case['plotdir'] = None  # if None, will not make plots
        #case['plotdir'] = plotdir
        case['setplot_file'] = 'setplot.py'

        caselist.append(case)

    return caselist
