import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def subplot_format(sp):
    xlim = [-3,3]
    ylim = [-3,3]

    plt.xlim(xlim)
    plt.ylim(ylim)

    sp.tick_params(
        axis='both',       # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        left='off',      # ticks along the bottom edge are off
        right='off',
        labelleft='off',# ticks along the top edge are off
        labelbottom='off') # labels along the bottom edge are off

def plot_laplacian(X_true, Y_true, lap_true, laplace_inf, xloc_syn, yloc_syn, skl_div_syn, T, t_train = [0], state = 'var'):
    assert len(T) < len(X_true)

    # Settings
    levels = 50
    if state == 'var':
        min_kl_cb = -180.
    else:
        min_kl_cb = np.min(np.array(skl_div_syn)[T])


    names = ["True Laplacian", "Inferred Laplacian", "Signed KL divergence"]

    # Settings: Plot
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    # Generate plot data
    levels_min = [np.min(lap_true), np.min(laplace_inf), min_kl_cb]
    levels_max = [np.max(lap_true), np.max(laplace_inf), np.max(skl_div_syn)]
    contour_x  = [X_true, xloc_syn, xloc_syn]
    contour_y  = [Y_true, yloc_syn, yloc_syn]
    contour_z  = [lap_true, laplace_inf, skl_div_syn]

    # Create plot
    for row in range(3):
        fig = plt.figure(figsize=(15,3))

        for col, t in enumerate(T):
            sp = plt.subplot(1,len(T),col+1)

            # Plot specific data
            levs_skl = np.linspace(levels_min[row], levels_max[row], levels)
            if row == 0:
                skl = plt.contourf(contour_x[row], contour_y[row], contour_z[row][t], levels = levs_skl, alpha=0.7)
            else:
                skl = plt.contourf(contour_x[row][t], contour_y[row][t], contour_z[row][t], levels = levs_skl, alpha=0.7)

            if state == 'stat':
                x = [1.5,-1.5]; y = [0.,0.]

            if state == 'var':
                # Include centre of min laplacian
                i = np.argmin(lap_true[t])
                x = np.ndarray.flatten(X_true)[i]
                y = np.ndarray.flatten(Y_true)[i]
                if x>0:
                    x = 1.5; y = 0.
                else:
                    x = -1.5; y=0

            if state == 'rot':
                mu1 = lambda x: np.array([[1.5*np.cos(x)], [1.5*np.sin(x)]])
                mu2 = lambda x: np.array([[-1.5*np.cos(x)], [1.5*np.sin(-x)]])
    
                time = t_train[0][t]
                A1 = mu1(time)
                A2 = mu2(time)
                x = np.hstack([A1,A2]).T[:,0]
                y = np.hstack([A1,A2]).T[:,1]
    
            # True point
            plt.scatter(x,y,color='w',s=200,alpha=1.0,marker='o')
            plt.scatter(x,y,color='k',s=50,alpha=1.0,marker='o')
            # Formatting
            subplot_format(sp)
            cax = fig.add_axes([1.01, 0.05, 0.02, 0.9])

            # Colour bar
            ticks = np.linspace(levels_min[row] + 0.1*(levels_max[row]-levels_min[row]), levels_max[row] - 0.1*(levels_max[row]-levels_min[row]),8)
            cb = fig.colorbar(skl, cax=cax,format="%1.2f",ticks=ticks)
            ticklabs = cb.ax.get_yticklabels()
            cb.ax.set_yticklabels(ticklabs,ha='right')
            cb.ax.yaxis.set_tick_params(pad=80)  # your number may vary
            cb.ax.tick_params(labelsize=20)

            # Y label
            if t == T[0]:
                sp.text(-4, 0, names[row], va='center', rotation='vertical', fontsize=19)
            if row == 0:
                sp.text(0., 3.4, 'Frame: %d' % t, ha='center', rotation='horizontal', fontsize=20)

        plt.tight_layout()
    plt.show()