



def print_corr(tica, feat, ic=0, num=20, save=False, outname=None):
    ''' helpful printing function for finding correlations between features and tICs

    .. math::
    \mathbf{Corr}(X_i - \mu_i, \mathbf{\theta}_j) = \frac{1}{\sigma_{X_i - \mu_i}}\sum_l \sigma_{(X_i - \mu_i)(X_l - \mu_l} \mathbf{U}_{li}
    
    INPUTs
    
    tica: a pyemma tica object
    feat: a pyemma featurizer object with a featurizer.describe() method
    ic: Independent component index to find correlations to
    num: number of top features to print
    save: save to file or just print. Boolean
    outname: name of output file

    OUTPUT
    Feature label
    correlation value from correlation matrix
    feature index
    '''
    top = np.argsort(np.absolute(tica.feature_TIC_correlation[:,ic]))[::-1]
    out = [(feat.describe()[top[i]], tica.feature_TIC_correlation[:,ic][top[i]], top[i]) for i in range(num)]
    if save==True:
       np.savetxt(f'{outname}.txt', out)
    print('Feature Label       Correlation      Feature index')
    for line in out: print(line)


def plot_by_index(ind, ylabel, savename): 
     plt.plot(time, cat_a[:, ind], label='chain A') 
     plt.plot(time, cat_b[:, ind], label='chain B') 
     plt.xlabel(f'Time (\u03BCs)', fontsize=16) 
     plt.ylabel(ylabel) 
     plt.legend() 
     plt.savefig(savename) 
     plt.clf() 

