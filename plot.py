import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def r_band():
    # R Band figure
    fr, ar = plt.subplots(2, 2, figsize=(10, 10), sharex=False, sharey=False)
    fr.suptitle('R Band Differences')
    #   Gini Coefficient
    rg = sns.distplot(ro[(ro.G > -0.95) & (ro.G < -0.5)]['G'].values, bins=100, color='r', ax=ar[0, 0],
                      hist=False, kde=True, label='type1')
    rg.set(xlabel='Gini Coefficient', ylabel='Gaussian Kernel Density')

    sns.distplot(rt[(rt.G > -0.95) & (rt.G < -0.5)]['G'].values, bins=100, color='b', ax=ar[0, 0],
                 hist=False, kde=True, label='type2')
    #   Concentration Index
    rc = sns.distplot(ro['C'].values, bins=100, color='m', ax=ar[0, 1],
                      hist=False, kde=True, label='type1')
    rc.set(xlabel='Concentration Index')
    sns.distplot(rt['C'].values, bins=100, color='y', ax=ar[0, 1],
                 hist=False, kde=True, label='type2')
    #   Moment Index
    rm = sns.distplot(ro[(ro.M > -3.5) & (ro.M < -1.2)]['M'].values, bins=100, color='g', ax=ar[1, 0],
                      hist=False, kde=True, label='type1')
    rm.set(xlabel='Moment Index', ylabel='Gaussian Kernel Density')
    sns.distplot(rt[(rt.M > -3.5) & (rt.M < -1.2)]['M'].values, bins=100, color='c', ax=ar[1, 0],
                 hist=False, kde=True, label='type2')
    #   Asymmetry Index
    rm = sns.distplot(ro[(ro.A > 0) & (ro.A < 0.8)]['A'].values, bins=100, color='brown', ax=ar[1, 1],
                      hist=False, kde=True, label='type1')
    rm.set(xlabel='Asymmetry Index')
    sns.distplot(rt[(rt.A > 0) & (rt.A < 0.8)]['A'].values, bins=100, color='hotpink', ax=ar[1, 1],
                 hist=False, kde=True, label='type2')


def g_band():
    # G Band figure
    fg, ag = plt.subplots(2, 2, figsize=(10, 10), sharex=False, sharey=False)
    fg.suptitle('G Band Differences')
    #   Gini Coefficient
    gg = sns.distplot(go[(go.G > -0.95) & (go.G < -0.5)]['G'].values, bins=100, color='r', ax=ag[0, 0],
                      hist=False, kde=True, label='type1')
    gg.set(xlabel='Gini Coefficient', ylabel='Gaussian Kernel Density')

    sns.distplot(gt[(gt.G > -0.95) & (gt.G < -0.5)]['G'].values, bins=100, color='b', ax=ag[0, 0],
                 hist=False, kde=True, label='type2')
    #   Concentration Index
    gc = sns.distplot(go['C'].values, bins=100, color='m', ax=ag[0, 1],
                      hist=False, kde=True, label='type1')
    gc.set(xlabel='Concentration Index', ylim=(0, 4.5))
    sns.distplot(gt['C'].values, bins=100, color='y', ax=ag[0, 1],
                 hist=False, kde=True, label='type2')
    #   Moment Index
    gm = sns.distplot(go[(go.M > -3.5) & (go.M < -1.2)]['M'].values, bins=100, color='g', ax=ag[1, 0],
                      hist=False, kde=True, label='type1')
    gm.set(xlabel='Moment Index', ylabel='Gaussian Kernel Density', ylim=(0, 2.5))
    sns.distplot(gt[(gt.M > -3.5) & (gt.M < -1.2)]['M'].values, bins=100, color='c', ax=ag[1, 0],
                 hist=False, kde=True, label='type2')
    #   Asymmetry Index
    gm = sns.distplot(go[(go.A > 0) & (go.A < 0.8)]['A'].values, bins=100, color='brown', ax=ag[1, 1],
                      hist=False, kde=True, label='type1')
    gm.set(xlabel='Asymmetry Index', ylim=(0, 6))
    sns.distplot(gt[(gt.A > 0) & (gt.A < 0.8)]['A'].values, bins=100, color='hotpink', ax=ag[1, 1],
                 hist=False, kde=True, label='type2')


# get the data
# r ---> 'r band',o/t --->type 1/2
ro = pd.read_csv('type1r.csv')
rt = pd.read_csv('type2r.csv')
go = pd.read_csv('type1g.csv')
gt = pd.read_csv('type2g.csv')
sns.set(style='darkgrid', palette='muted', color_codes=True, font_scale=1.5)
r_band()
g_band()
# sns.jointplot(ro['C'], rt['C'], kind='kde', color='b')
plt.show()
