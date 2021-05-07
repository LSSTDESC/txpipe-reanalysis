import numpy as np
import pylab as mplot
import os
import scipy
import pandas as pd
from getdist import plots, MCSamples
from chainconsumer import ChainConsumer
import matplotlib.pyplot as plt

import sys

sys.path.append('/Library/Python/2.7/site-packages')

from chainconsumer import ChainConsumer

Color = ['#d45e00', 'grey','k', 'purple']

font = {'size'   : 18}
mplot.rc('font', **font)
mplot.rc('text', usetex=True)
mplot.rc('font', family='serif')

def load_cosmosis_chain(filename,headers,caster=None,doS8=False,flipz=False,add_cols_end=['prior','like','posterior','weights'],add_cols_front=None,rename_cols=None):
    ''' Load a cosmosis chain.

    Parameters
    ==========
    filename : str,         
        cosmosis chain filename
    headers : list,
        list of the headers in the file
    caster : str, default = None
        optionally cast a column with name given to a float from a str
    doS8 : bool, default = False
        optionally compute S8 as sigma8*(omega_m/0.3)**0.5
    flipz : bool, default = False 
        optionally flip the sign on any redshift bias values
    add_cols_end: list, default = ['prior','like','posterior',weight']
        optionally add columns to the end of the header
    add_cols_front: list, default = None
        optinally add columns to the front of the header
    rename_cols: dictionary, default = None
        optionally rename the column names, if passing give a mapper object to pass to pandas 
        

    '''
    if add_cols_end:
        cosmosis_headers = headers + add_cols_end
    if add_cols_front:
        cosmosis_headers = add_cols_front + headers
    cosmosis_data = pd.read_csv(filename,delim_whitespace=True,comment='#',names=cosmosis_headers)
    if caster:
        cosmosis_data[caster] = cosmosis_data[caster].astype(float)
    if flipz:
        for colname in cosmosis_headers:
            if 'z' in colname:
                cosmosis_data[colname] = -cosmosis_data[colname]
    if rename_cols:
        cosmosis_data = cosmosis_data.rename(columns=rename_cols)
    if doS8:
        cosmosis_data['S8'] = cosmosis_data['sigma8']*(cosmosis_data['Omega_m']/0.3)**0.5
    return cosmosis_data
        
        
def load_samples(chains,ids,parameters):
    ''' Load a cosmosis chain into a list of getdist MCSamples objects.

    Parameters
    ==========
    chains : list of cosmosis data,         
        cosmosis chain filename
    ids : list,
        list of the ID's of the chains, e.g. ['desy1','hsc']
    parameters : list, 
        list of the parameters to include in the samples, e.g. ['omega_m','sigma8']
    
    Returns
    ==========
    Samples : list, 
        list of MCSample objects
        
    '''
    
    if type(chains)==pd.core.frame.DataFrame:
        data = chains
        paramdata = data.drop(columns=['prior','like','posterior','weights'],errors='ignore')
        drop_columns = [col for col in paramdata.columns if col not in parameters]
        paramdata = paramdata.drop(columns=drop_columns)
        samples = MCSamples(samples=paramdata.values,weights=data['weights'].values,
                    names = paramdata.columns,
                    labels = paramdata.columns)
        return samples
    else:
        Samples = []
        for i in range(len(chains)):
            data = chains[i]
            paramdata = data.drop(columns=['prior','like','posterior','weights'],errors='ignore')
            drop_columns = [col for col in paramdata.columns if col not in parameters]
            paramdata = paramdata.drop(columns=drop_columns)
            samples = MCSamples(samples=paramdata.values,weights=data['weights'].values,
                        names = paramdata.columns,
                        labels = paramdata.columns)

            Samples.append(samples)
    return Samples

def triangle_plot(samples,ids,parameters,width=5.0):
    ''' Plot a getdist triangle plot.

    Parameters
    ==========
    samples : list,         
        list of MCSamples
    ids : list, 
        list of the ID's of the samples, e.g. ['desy1','hsc']
    parameters : list,
        list of the parameters to plot, e.g. ['omega_m','sigma8']
    
    Returns
    ==========
    Samples : list, 
        list of MCSample objects
        
    '''
    g = plots.get_subplot_plotter(width)
    g.settings.figure_legend_frame = False
    g.settings.alpha_filled_add=0.4
    g.settings.title_limit_fontsize = 14
    g.triangle_plot(samples, parameters, 
        filled=True, 
        legend_labels=ids, 
        legend_loc='upper right', 
        line_args=[{'ls':'--', 'color':'green'},
                   {'lw':2, 'color':'darkblue'}], 
        contour_colors=['green','darkblue'],
        title_limit=1, # first title limit (for 1D plots) is 68% by default
        markers={'x2':0})

class ParLims():
    """ A class to store parameter values and limits.
    Work in progress, currently just takes a chainconsumer summary, prints it out. Would like to add 
    Latex, table capabilities, etc. """ 
    
    def __init__(self, summary):
        self.summary = summary 
        
    def show_pm(self):
        for parameter in self.summary.keys():
            upper = self.summary[parameter][2]-self.summary[parameter][1]
            lower = self.summary[parameter][1]-self.summary[parameter][0]
            mid = self.summary[parameter][1]

            print('{parameter} : {mid:.3f}, + {upper:.3f} - {lower:.3f}'.format(parameter=parameter,upper=upper,lower=lower,mid=mid))
    
    def get_pm(self,parameter):
        upper = self.summary[parameter][2]-self.summary[parameter][1]
        lower = self.summary[parameter][1]-self.summary[parameter][0]
        mid = self.summary[parameter][1]
        
        return mid, lower, upper
    
def get_parlims_mean_getDist(chain,parameters,interval='68'):
    """ Get the parameter values and 68% upper/lower bounds for input parameters at the mean of the posterior.  
    
    Parameters
    ==========
    chain : pandas.dataframe,         
        pandas dataframe holding a cosmosis chain
    parameters : list, 
        list of the parameters to get statistics for, e.g. ['omega_m','sigma8']
    statistic : str,
        type of parameter, e.g. 'mean' to apply to the posterior to determine the central parameter value
    
    Returns
    ==========
    mean : float, 
        mean parameter value
    upper_limit : float,
        upper parameter limit
    lower_limit : float,
        lower parameter limit
    """
    samples = load_samples(chain,['chain'],parameters)
    
    res = {}
    interval = '68'
    if interval == '68':
        index = 0
    elif interval == '95':
        index = 1
    elif interval =='99':
        index = 3
    else:
        RaiseExcept('Please specify a 68, 95, or 99 limit')
    for parameter in parameters:
        lim_stats = samples.getMargeStats()
        lim_par_stats = lim_stats.parWithName(parameter)
        mean = lim_par_stats.mean
        lims = lim_par_stats.limits[index] # get the 68
        lower, upper = lims.lower, lims.upper
        res[parameter] = [lower, mean, upper]
    
    parlims = ParLims(res)
    return parlims


def get_parlims_maxpost_chainconsumer(chain,parameters,kde=False):
    """ Get the parameter values and 68% upper/lower bounds for input parameters at the maximum of the posterior.
    
    Parameters
    ==========
    chain : pandas.dataframe,         
        pandas dataframe holding a cosmosis chain
    parameters : list, 
        list of the parameters to get statistics for, e.g. ['omega_m','sigma8']
    statistic : str,
        type of parameter, e.g. 'mean' to apply to the posterior to determine the central parameter value
    
    Returns
    ==========
    mean : float, 
        mean parameter value
    upper_limit : float,
        upper parameter limit
    lower_limit : float,
        lower parameter limit
        
    """
    c = ChainConsumer()
    w = chain.weights
    drop_columns = [col for col in chain.columns if col not in parameters]
    chain = chain.drop(columns=drop_columns)
    
    c.add_chain(chain.values,parameters=list(chain.columns),weights=w)

    c.configure(statistics='max',kde=kde,sigma2d=False)
    
    summary = c.analysis.get_summary() 
    
    parlims = ParLims(summary)
    
    
    return parlims
    

def plot_parlims(parlims,values,parameters,title_parameters,shade=False,parlim_labels=None,stats=[]):
    if len(stats)==1:
        parlims = [parlims]
    deljs = []
    colors = ['tab:blue','tab:orange','tab:green']
    fmts = ['o','d','X']
    alphas = [1.0,0.5]
    k = 0
    ys = []
    fig, axes = plt.subplots(1, len(parameters)+1,sharey=True)
    for parlim in parlims:
        i = 0
        j = 0
        for par in parameters:
            i = 0
            for survey_parlims in parlim:
                mid,lower,upper = survey_parlims.get_pm(par)
                err = np.array([[lower,upper]]).reshape(2,1)
                axes[j].errorbar([mid],[i+k*0.1],xerr=err,fmt=fmts[k],alpha=alphas[k],elinewidth=3.0,color=colors[j])
                ys.append(i)
                i+=1
                if shade==True:
                    axes[j].axvspan(mid-lower, mid+upper, alpha=0.5, color=colors[i])
            axes[j].title.set_text(title_parameters[j]) 
            j = j+1
        axes[j].errorbar([mid],[i+k*0.1],xerr=err,fmt=fmts[k],alpha=alphas[k],elinewidth=3.0,color='black',label=stats[k])
        plt.legend(loc='upper right' , bbox_to_anchor=(1.2, 1.2),  )
        k = k+1
        deljs.append(j)
    #for jj in deljs:
    #    axes[jj].set_visible(False)
    plt.yticks(range(len(parameters)),values)
    plt.show()
    