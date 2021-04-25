def many_plot(alphas = [0.3, 0.6, 0.9], stim_a = stimulusA , stim_x = stimulusX, y = unconditioned):
    """Plotting function for getting the five rows of the final plot. 

    Args:
        alphas (list, optional): List of learning rates. Defaults to [0.3, 0.6, 0.9].
        stim_a (np.array, optional): Numpy array of stimulus "A". Defaults to stimulusA.
        stim_x (np.array, optional): Numpy array of stimulus "X". Defaults to stimulusX.
        y (np.array, optional): Numpy array of unconditioned stimuli. Defaults to unconditioned.

    Returns:
        fig: Returns the plot.
    """    
    fig, ax = plt.subplots(5,1, figsize = (12,24))
    ax[0].plot(stim_a, color = "C7")
    ax[1].plot(stim_x, color = "C8")
    ax[2].plot(y, color = "C9")
    ax[0].set(xlabel = "time step", ylabel = "state", title = "Stimulus A")
    ax[1].set(xlabel = "time step", ylabel = "state", title = "Stimulus X")
    ax[2].set(xlabel = "time step", ylabel = "state", title = "Unconditioned Stimulus")

    for alpha in alphas:
        placeholder = MostWagner(y, stim_a, stim_x, alpha)
        placeholder.run_sim()
        a,x,_ = list(zip(*placeholder.data))
        ax[3].plot(a, label = f"{alpha}")
        ax[4].plot(x, label = f"{alpha}")
    ax[3].legend(loc = "upper left")
    ax[4].legend(loc = "upper left")
    ax[3].set(xlabel = "time step", ylabel = "asso.", title = r"$V_A$")
    ax[4].set(xlabel = "time step", ylabel = "asso.", title = r"$V_X$")
    ax[4].set_ylim(0,1)
    ax[3].set_ylim(0,1)

    fig.tight_layout(pad=2.0)
    
    return fig
