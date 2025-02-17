
def setplot(plotdata, case={}):

    """
    Run the code for making fgmax and gauge plots and fgout animations.
    This version does not return plotdata for making time frame plots.
    Add code to set plotdata and return it if desired.
    """

    location = 'Newport' # need to change module names below for other locations

    outdir = case['outdir']
    plotdir = case['plotdir']
    dtopofiles = case['dtopofiles']

    dtopofile = dtopofiles[0][-1]
    event = os.path.splitext(os.path.split(dtopofile)[-1])[0]

    run_name = '%s_%s' % (location,event)

    if 1:
        gauges_plotdir = plotdir + '/gauges'
        os.system('mkdir -p %s' % gauges_plotdir)
        plot_gauges_site.make_all_plots_and_report(outdir, gauges_plotdir,
                         location=location, event=event,
                         gaugenos='all', sea_level=0.)

    if 1:
        try:
            fgmax_plotdir = plotdir + '/fgmax'
            os.system('mkdir -p %s' % fgmax_plotdir)
            process_fgmax.make_all_fgmax_plots(outdir, fgmax_plotdir,
                                           location=location, event=event)
        except:
            print('*** problem with process_fgmax in %s, skipping' \
                  % run_name)
    if 1:
        make_fgout_animation.make_anim(outdir, plotdir, location, event)

    if 0:
        make_html_index(plotdir,event)

    plotdata = None

    # add code to set plotdata (as in other setplot.py modules) if you also
    # want to make time frame plots.

    return plotdata
