## ---------------------------
##
##  Load the custom BU code for covasim
##
## Authors: Brian Gregor, Wenrui Li 
##          Boston University
##
## Date Created: 2020-07-31
##
## Email: bgregor@bu.edu
##
## ---------------------------
# This MUST COME FIRST in order to monkeypatch the covasim People and MultiSim classes.

import covasim as cv
import covasim.run as cr





# =============================================================================
# Tweak the People class to fit our needs
# =============================================================================
# Add some people characteristics to the PeopleMeta class
BU_attrs = ['campus','undergrad','test_cat','category', 'campResident', 'full_info_id','labTestEveryNHours']  
cv.defaults.PeopleMeta.person += BU_attrs
cv.defaults.PeopleMeta.all_states += BU_attrs

# =============================================================================
# Update the Multisim plot() function to properly pass on the 
# show_args argument to the main plotting function.
# =============================================================================

def new_plot(self, to_plot=None, inds=None, plot_sims=False, color_by_sim=None, max_sims=5, colors=None, labels=None, alpha_range=None, plot_args=None, show_args=None, **kwargs):
    # Plot a single curve, possibly with a range
    if not plot_sims and self.which in ['combined', 'reduced']:
        fig = self.base_sim.plot(to_plot=to_plot, colors=colors, show_args=show_args, **kwargs)

    # PLot individual sims on top of each other
    else:

        # Initialize
        fig = None
        orig_show    = kwargs.get('do_show', True)
        orig_setylim = kwargs.get('setylim', True)
        kwargs['legend_args'] = sc.mergedicts({'show_legend':True}, kwargs.get('legend_args')) # Only plot the legend the first time

        # Handle indices
        if inds is None:
            inds = np.arange(len(self.sims))
        n_sims = len(inds)

        # Handle what style of plotting to use:
        if color_by_sim is None:
            if n_sims <= max_sims:
                color_by_sim = True
            else:
                color_by_sim = False

        # Handle what to plot
        if to_plot is None:
            if color_by_sim:
                to_plot = cvd.get_scen_plots()
            else:
                to_plot = cvd.get_sim_plots()

        # Handle colors
        if colors is None:
            if color_by_sim:
                colors = sc.gridcolors(ncolors=n_sims)
            else:
                colors = [None]*n_sims # So we can iterate over it
        else:
            colors = [colors]*n_sims # Again, for iteration

        # Handle alpha if not using colors
        if alpha_range is None:
            if color_by_sim:
                alpha_range = [0.8, 0.8] # We're using color to distinguish sims, so don't need alpha
            else:
                alpha_range = [0.8, 0.3] # We're using alpha to distinguish sims
        alphas = np.linspace(alpha_range[0], alpha_range[1], n_sims)

        # Plot
        for s,ind in enumerate(inds):
            sim = self.sims[ind]

            final_plot = (s == n_sims-1) # Check if this is the final plot

            # Handle the legend and labels
            if final_plot:
                merged_show_args  = show_args
                kwargs['do_show'] = orig_show
                kwargs['setylim'] = orig_setylim
            else:
                merged_show_args  = False # Only show things like data the last time it's plotting
                kwargs['do_show'] = False # On top of that, don't show the plot at all unless it's the last time
                kwargs['setylim'] = False

            # Optionally set the label for the first max_sims sims
            if labels is None and color_by_sim is True and s<max_sims:
                merged_labels = sim.label
            elif final_plot:
                merged_labels = labels
            else:
                merged_labels = ''

            # Actually plot
            merged_plot_args = sc.mergedicts({'alpha':alphas[s]}, plot_args) # Need a new variable to avoid overwriting
            fig = sim.plot(fig=fig, to_plot=to_plot, colors=colors[s], labels=merged_labels, plot_args=merged_plot_args, show_args=merged_show_args, **kwargs)

    return fig


cr.MultiSim.plot = new_plot

from .misc  import *
from .housing_and_classrooms import * 
from .parallel  import *
from .snapshots import *
from .interventions import *
from .exception import *
