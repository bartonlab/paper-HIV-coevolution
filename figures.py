"""
This Python file contains code for generating plots in:
    Parallel HIV-1 evolutionary dynamics in humans and rhesus macaques who develop broadly neutralizing antibodies
        by Kai Shimagaki
           Rebecca Lynch
           John Barton (jpbarton@pitt.edu)

Additional code to pre-process the data and pass it to these plotting
routines is contained in the Jupyter notebook `figures.ipynb`.
"""

#############  PACKAGES  #############

import sys, os
from copy import deepcopy

import numpy as np

import scipy as sp
import scipy.stats as st
import scipy.optimize as so

import pandas as pd

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg

from PIL import Image 

import seaborn as sns

from colorsys import hls_to_rgb

import mplot as mp


############# PARAMETERS #############

# GLOBAL VARIABLES

NUC = ['-', 'A', 'C', 'G', 'T']
REF = NUC[0]
EXT = '.pdf'

## Code Ocean directories
#MLB_DIR = '../data/HIV'
#WFS_DIR = '../data/wfsim'
#HIV_MPL_DIR = '../data/HIV/MPL'
#SIM_MPL_DIR = '../data/simulation/MPL'
#HIV_DIR = '../data/HIV'
#SIM_DIR = '../data/simulation'
#FIG_DIR = '../results'

# GitHub directories
MLB_DIR = 'src/Matlab'
WFS_DIR = 'src/wfsim'
HIV_MPL_DIR = 'src/MPL/HIV'
SIM_MPL_DIR = 'src/MPL/out'
HIV_DIR = 'data/HIV'
SIM_DIR = 'data/simulation'
FIG_DIR = 'figures/'

# Standard color scheme

BKCOLOR  = '#252525'
LCOLOR   = '#969696'
C_BEN    = '#EB4025' #'#F16913'
C_BEN_LT = '#F08F78' #'#fdd0a2'
C_NEU    =  LCOLOR   #'#E8E8E8' # LCOLOR
C_NEU_LT = '#E8E8E8' #'#F0F0F0' #'#d9d9d9'
C_DEL    = '#3E8DCF' #'#604A7B'
C_DEL_LT = '#78B4E7' #'#dadaeb'
C_MPL    = '#FFB511'
C_SL     = C_NEU
C_BNAB   = '#EA36F7'

# Plot conventions

def cm2inch(x): return float(x)/2.54
SINGLE_COLUMN   = cm2inch(8.8)
ONE_FIVE_COLUMN = cm2inch(11.4)
DOUBLE_COLUMN   = cm2inch(18.0)

GOLDR        = (1.0 + np.sqrt(5)) / 2.0
TICKLENGTH   = 3
TICKPAD      = 3
AXWIDTH      = 0.4

# paper style
FONTFAMILY   = 'Arial'
SIZESUBLABEL = 8
SIZELABEL    = 6
SIZETICK     = 6
SMALLSIZEDOT = 6.
SIZELINE     = 0.6

# grant style
#FONTFAMILY   = 'Arial'
#SIZESUBLABEL = 8
#SIZELABEL    = 8
#SIZETICK     = 8
#SMALLSIZEDOT = 6. * 1.3
#SIZELINE     = 0.6

# slides style
#FONTFAMILY   = 'Avenir'
#SIZESUBLABEL = 20
#SIZELABEL    = 20
#SIZETICK     = 20
#SMALLSIZEDOT = 6. * 7
#SIZELINE     = 1.5
#SLIDE_WIDTH  = 10.5


FIGPROPS = {
    'transparent' : True,
    #'bbox_inches' : 'tight'
}

DEF_BARPROPS = {
    'lw'          : SIZELINE/2,
    'width'       : 0.25,
    'edgecolor'   : BKCOLOR,
    'align'       : 'center', #other option: edge
    'orientation' : 'vertical'
}

DEF_HISTPROPS = {
    'histtype'    : 'bar',
    'lw'          : SIZELINE/2,
    'rwidth'      : 0.8,
    'ls'          : 'solid',
    'edgecolor'   : BKCOLOR,
    'alpha'       : 0.5
}

DEF_ERRORPROPS = {
    'mew'        : AXWIDTH,
    'markersize' : SMALLSIZEDOT/2,
    'fmt'        : 'o',
    'elinewidth' : SIZELINE/2,
    'capthick'   : 0,
    'capsize'    : 0
}

DEF_LINEPROPS = {
    'lw' : SIZELINE,
    'ls' : '-'
}

DEF_LABELPROPS = {
    'family' : FONTFAMILY,
    'size'   : SIZELABEL,
    'color'  : BKCOLOR
}

DEF_SUBLABELPROPS = {
    'family'  : FONTFAMILY,
    'size'    : SIZESUBLABEL+1,
    'weight'  : 'bold',
    'ha'      : 'center',
    'va'      : 'center',
    'color'   : 'k',
    'clip_on' : False
}

DEF_TICKLABELPROPS = {
    'family' : FONTFAMILY,
    'size'   : SIZETICK,
    'color'  : BKCOLOR
}

DEF_TICKPROPS = {
    'length'    : TICKLENGTH,
    'width'     : AXWIDTH/2,
    'pad'       : TICKPAD,
    'axis'      : 'both',
    'direction' : 'out',
    'colors'    : BKCOLOR,
    'bottom'    : True,
    'left'      : True,
    'top'       : False,
    'right'     : False
}

DEF_MINORTICKPROPS = {
    'length'    : TICKLENGTH-1.25,
    'width'     : AXWIDTH/2,
    'axis'      : 'both',
    'direction' : 'out',
    'which'     : 'minor',
    'color'     : BKCOLOR
}

DEF_AXPROPS = {
    'linewidth' : AXWIDTH,
    'linestyle' : '-',
    'color'     : BKCOLOR
}

PARAMS = {'text.usetex': False, 'mathtext.fontset': 'stixsans', 'mathtext.default' : 'regular'}
plt.rcParams.update(PARAMS)


############# PLOTTING  FUNCTIONS #############


def plot_figure_ch505_ch848_circle():
    """
    CH505 and CH848 inferred selection coefficients.
    """

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN # SLIDE_WIDTH
    hshrink = 0.50
    goldh   = w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    box_top = 0.93
    box_bot = 0.03
    dy      = 0 #0.08  
    y_circ  = (box_top - box_bot - dy) * (5/5)  # decide the relative height of the circle plots
    y_enr   = (box_top - box_bot - dy) * (0/5)  # decide the relative height of the enrichment plots
    x_circ  = y_circ * goldh / w  # ensure axes are square
    dx      = (1 - 2*x_circ) / 3  # position circle plots in center
    x0      = -0.02

    box_circ_505 = dict(left=x0 + 1*dx + 0*x_circ, right=x0 + 1*dx + 1*x_circ, bottom=box_top-y_circ, top=box_top)
    box_circ_848 = dict(left=x0 + 2*dx + 1*x_circ, right=x0 + 2*dx + 2*x_circ, bottom=box_top-y_circ, top=box_top)
    # box_enr_505 = dict(left=1*dx + 0*x_circ, right=1*dx + 1*x_circ, bottom=box_top-y_circ-y_enr-dy, top=box_top-y_circ-dy)
    # box_enr_848 = dict(left=2*dx + 1*x_circ, right=2*dx + 2*x_circ, bottom=box_top-y_circ-y_enr-dy, top=box_top-y_circ-dy)

    gs_circ_505 = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_circ_505)
    gs_circ_848 = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_circ_848)
    # gs_enr_505 = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_enr_505)
    # gs_enr_848 = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_enr_848)

    ax_circ_505 = plt.subplot(gs_circ_505[0, 0])
    ax_circ_848 = plt.subplot(gs_circ_848[0, 0])
    # ax_enr_505 = plt.subplot(gs_enr_505[0, 0])
    # ax_enr_848 = plt.subplot(gs_enr_848[0, 0])
    
    #dx = -0.08
    dy_sub   = 0.02
    dy_label = 0.035
    dx_label = 0.033

    ## a -- circle plot

#    # _ SLIDES
#    w       = 8
#    fig     = plt.figure(figsize=(w, w))
#    ax_circ = plt.subplot(111)
#    # ^ SLIDES

    label2ddr = { 'tat exon 2': 0.13,
                  'rev exon 2': 0.10 }

    ch505_s     = 'data/703010505-3-poly.csv'
    ch505_index = 'data/703010505-3-index.csv'

    sig_s, sig_site_real, sig_nuc_idx = plot_circle(ax_circ_505, ch505_s, ch505_index, label2ddr)

    ax_circ_505.text(box_circ_505['left']+dy_sub, box_circ_505['top']+dy_sub, 'a'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_circ_505.text(box_circ_505['left']+dx_label, box_circ_505['top']-dy_label, 'CH505', transform=fig.transFigure, **DEF_LABELPROPS)


    label2ddr = { 'tat exon 2': 0.10,
                  'rev exon 2': 0.15 }

    ch848_s     = 'data/703010848-3-poly.csv'
    ch848_index = 'data/703010848-3-index.csv'

    sig_s, sig_site_real, sig_nuc_idx = plot_circle(ax_circ_848, ch848_s, ch848_index, label2ddr)

    ax_circ_848.text(box_circ_848['left']+dy_sub, box_circ_848['top']+dy_sub, 'b'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_circ_848.text(box_circ_848['left']+dx_label, box_circ_848['top']-dy_label, 'CH848', transform=fig.transFigure, **DEF_LABELPROPS)

#    # _ SLIDES
#    plt.savefig('%s/new-slides-cap256-circle.pdf' % FIG_DIR, dpi=1000, **FIGPROPS)
#    plt.savefig('%s/new-slides-cap256-circle.png' % FIG_DIR, dpi=1000, **FIGPROPS)
#    plt.close(fig)
#    w       = 8
#    fig     = plt.figure(figsize=(w, 2.85))
#    ax_smpl = plt.subplot(311)
#    # ^ SLIDES

#     dx =  0.04
#     dy = -0.02
#     ax_circ.text(0.06, box_circ['top']+dy, 'b'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

#     ## a -- trajectory plot

#     pprops = { 'xticks':      [30, 190, 350, 510, 670],
#                'yticks':      [0, 1],
#                'yminorticks': [0.25, 0.5, 0.75],
#                'nudgey':      1.1,
#                'xlabel':      'Time (days)',
#                'ylabel':      'Variant frequency\nin VRC26 epitope',
#                'plotprops':   {'lw': 1.5*SIZELINE, 'ls': '-', 'alpha': 0.75 },
#                'axoffset':    0.1,
#                'theme':       'open' }

#     for i in range(len(var_tag)-1):
#         xdat = [times]
#         ydat = [var_traj[i]]
#         mp.line(ax=ax_traj, x=xdat, y=ydat, colors=[var_c[i]], **pprops)

#     xdat = [times]
#     ydat = [var_traj[len(var_tag)-1]]
#     mp.plot(type='line', ax=ax_traj, x=xdat, y=ydat, colors=[var_c[len(var_tag)-1]], **pprops)

#     ax_traj.text(0.06, box_traj['top']+dy, 'a'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

#     invt = ax_traj.transData.inverted()
#     xy1  = invt.transform((  0, 0))
#     xy2  = invt.transform((7.5, 9))
#     xy3  = invt.transform((3.0, 9))

#     legend_dx1 = 1*(xy1[0]-xy2[0])
#     legend_dx2 = 1*(xy1[0]-xy3[0])
#     legend_dy  = 1*(xy1[1]-xy2[1])

#     traj_legend_x =  80
#     traj_legend_y =  0.95
#     traj_legend_d = -0.05
#     traj_legend_t = ['Superinfecting\nvariant', 'Other']
#     traj_legend_c = [C_MPL, C_NEU]
#     for k in range(len(traj_legend_t)):
#         mp.line(ax=ax_traj, x=[[traj_legend_x + legend_dx1, traj_legend_x + legend_dx2]],
#                 y=[[traj_legend_y + (1.5 * k * legend_dy), traj_legend_y + (1.5 * k * legend_dy)]],
#                 colors=[traj_legend_c[k]], plotprops=dict(lw=1.5*SIZELINE, ls='-', clip_on=False))
#         ax_traj.text(traj_legend_x, traj_legend_y + (1.5 * k * legend_dy), traj_legend_t[k], ha='left', va='center', **DEF_LABELPROPS)

#     ## c -- selection in the VRC26 epitope

#     site_rec_props  = dict(height=1, width=1, ec=None, lw=AXWIDTH/2, clip_on=False)
#     codon_rec_props = dict(height=4, width=3, ec=BKCOLOR, fc='none', lw=AXWIDTH/2, clip_on=False)
#     cBG             = '#F5F5F5'
#     rec_patches     = []
#     TF_dots_x       = [ 0.5]
#     TF_dots_y       = [-3.5]

#     sig_s         = np.array(sig_s)
#     sig_site_real = np.array(sig_site_real)
#     sig_nuc_idx   = np.array(sig_nuc_idx)
#     eidx          = 0
#     sub_box       = 0
#     sub_dx        = 3 + 0.3
#     for i in range(epitope_start[eidx], epitope_end[eidx]+1, 3):
#         for sub_i in range(3):
#             TF_dots_x.append(sub_box*sub_dx + sub_i + 0.5)
#             TF_dots_y.append(4-NUC.index(df_index.iloc[i+sub_i].TF)+0.5)
#             idxs     = sig_site_real==i+sub_i
#             temp_s   = sig_s[idxs]
#             temp_nuc = sig_nuc_idx[idxs]
#             for j in range(len(NUC)-1):
#                 if df_index.iloc[i+sub_i].TF==NUC[j+1]:
#                     continue
#                 c = cBG
#                 if j in temp_nuc:
#                     t = temp_s[list(temp_nuc).index(j)] / 0.05
#                     if np.fabs(t)>1:
#                         t /= np.fabs(t)
#                     if t>0:
#                         c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
#                     else:
#                         c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
#                 rec = matplotlib.patches.Rectangle(xy=(sub_box*sub_dx + sub_i, 3-j), fc=c, **site_rec_props)
#                 rec_patches.append(rec)
#         rec = matplotlib.patches.Rectangle(xy=(sub_box*sub_dx, 0), **codon_rec_props)
#         rec_patches.append(rec)
#         sub_box += 1

#     for i in range(-5, 5+1, 1):
#         c = cBG
#         t = i/5
#         if t>0:
#             c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
#         else:
#             c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
#         rec = matplotlib.patches.Rectangle(xy=(sub_box*sub_dx - 6.5 + i, -4), fc=c, **site_rec_props)
#         rec_patches.append(rec)

#     invt = ax_smpl.transData.inverted()
#     xy1  = invt.transform((0,0))
#     xy2  = invt.transform((0,9))
#     coef_legend_dy = (xy1[1]-xy2[1]) # multiply by 3 for slides/poster
#     c   = cBG
#     rec = matplotlib.patches.Rectangle(xy=(0, -4 + 4*coef_legend_dy), fc=c, **site_rec_props) # paper
# #    rec = matplotlib.patches.Rectangle(xy=(0, -4 + 6*coef_legend_dy), fc=c, **site_rec_props) # slides
#     rec_patches.append(rec)

#     for patch in rec_patches:
#         ax_smpl.add_artist(patch)

#     pprops = { 'colors': [BKCOLOR],
#                'xlim': [0, epitope_end[eidx]-epitope_start[eidx] + (sub_box * (sub_dx - 3)) + 1],
#                'ylim': [0, 4.1],
#                'xticks': [],
#                'yticks': [],
#                'plotprops': dict(lw=0, s=0.2*SMALLSIZEDOT, marker='o', clip_on=False),
#                'ylabel': '',
#                'theme': 'open',
#                'hide' : ['top', 'bottom', 'left', 'right'] }

#     mp.plot(type='scatter', ax=ax_smpl, x=[TF_dots_x], y=[TF_dots_y], **pprops)

#     ax_smpl.text(0.06, box_smpl['top']+0.02, 'c'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

#     txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
#     for i in range(len(NUC)-1):
#         ax_smpl.text(-0.85, 3-i+0.5, NUC[i+1], clip_on=False, **txtprops)

#     txtprops['ha'] = 'left'
#     ax_smpl.text(1.3, -3.5, 'TF nucleotide', clip_on=False, **txtprops)
#     ax_smpl.text(1.3, -3.5 + 4*coef_legend_dy, 'Not observed', clip_on=False, **txtprops) # paper
# #    ax_smpl.text(1.3, -3.5 + 6*coef_legend_dy, 'Not observed', clip_on=False, **txtprops) # slides

#     txtprops['ha'] = 'center'
#     txtprops['va'] = 'top'
#     for i in range((epitope_end[eidx]-epitope_start[eidx]+1)//3):
#         ax_smpl.text(1.5+i*sub_dx, -0.5, 160+i, clip_on=False, **txtprops)

#     ax_smpl.text(sub_box*sub_dx - 11, -4.5, -5, clip_on=False, **txtprops)
#     ax_smpl.text(sub_box*sub_dx -  6, -4.5,  0, clip_on=False, **txtprops)
#     ax_smpl.text(sub_box*sub_dx -  1, -4.5,  5, clip_on=False, **txtprops)
#     ax_smpl.text(sub_box*sub_dx -  5.5, -6.0, 'Inferred selection\ncoefficient, $\hat{s}$ (%)', clip_on=False, **txtprops)

# #    # _ SLIDES
# #    plt.savefig('%s/new-slides-cap256-selection.pdf' % FIG_DIR, dpi=1000, **FIGPROPS)
# #    plt.savefig('%s/new-slides-cap256-selection.png' % FIG_DIR, dpi=1000, **FIGPROPS)
# #    plt.close(fig)
# #    # ^ SLIDES

    # MAKE LEGEND

    invt = ax_circ_505.transData.inverted()
    xy1  = invt.transform((0,0))
    xy2  = invt.transform((0,15))

    coef_legend_x  = 0.90
    coef_legend_dx = 0.05
    coef_legend_y  = -0.85
    coef_legend_dy = xy1[1]-xy2[1]

    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    s_mult   = 40*SMALLSIZEDOT
    ex_s     = [-0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.10]
    show_s   = [   1,     0,     0,     0,     0, 1,    0,    0,    0,    0,    1]
    c_s      = [C_DEL, C_BEN]
    for i in range(len(ex_s)):
        plotprops = dict(lw=0, marker='o', s=np.fabs(ex_s[i])*s_mult, clip_on=False)
        mp.scatter(ax=ax_circ_505, x=[[coef_legend_x + i*coef_legend_dx]], y=[[coef_legend_y]], colors=[c_s[ex_s[i]>0]], plotprops=plotprops)
        if show_s[i]:
            ax_circ_505.text(coef_legend_x + i*coef_legend_dx, coef_legend_y + 0.75*coef_legend_dy, '%d' % (100*ex_s[i]), **txtprops)
    ax_circ_505.text(coef_legend_x + 5*coef_legend_dx, coef_legend_y + (2.25*coef_legend_dy), 'Inferred selection\ncoefficient, $\hat{s}$ (%)', **txtprops)

    # SAVE FIGURE

    plt.savefig('%sfig-ch505-ch848-circle%s' % (FIG_DIR, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('CH505-CH848 done.')


def plot_circle(ax, s_file, index_file, label2ddr):
    """
    Create a circle plot showing selection, large elements of the inverse covariance, and HIV sequence features.
    """

    df = pd.read_csv('%s' % (s_file), comment='#', memory_map=True)
    df_index = pd.read_csv('%s' % (index_file), comment='#', memory_map=True)

    # cov  = np.loadtxt('%s/covariance-%s-poly-seq2state.dat' % (HIV_MPL_DIR, tag))
    # num  = np.loadtxt('%s/numerator-%s-poly-seq2state.dat' % (HIV_MPL_DIR, tag))
    # cinv = np.linalg.inv(cov)
    # ds   = cinv / 1e4

    # Check for selection coefficients that are significant: |s| > mult * sigma_s

    s_min = 0.001
    n_sig = 0

    sig_s         = []
    sig_site_real = []
    sig_nuc_idx   = []

    #print('%s\nvariant\ts\t(sigma)' % (tag))
    for df_iter, df_entry in df.iterrows():
        if df_entry.nucleotide!='-':
            s_val = df_entry.s_MPL
            if np.fabs(s_val)>s_min and df_entry.nonsynonymous>0:
                n_sig += 1
                sig_s.append(s_val)
                sig_site_real.append(df_entry.alignment_index)
                sig_nuc_idx.append(NUC.index(df_entry.nucleotide) - 1)

#    print('')
#    print('%d/%d (%d%%) significant at %.2f sigma' % (n_sig, len(df_info), 100*n_sig/len(df_info), mult))
#    print('')

    # Sequence landmarks

    seq_range = [[ 790, 2292], [2085, 5096], [5041, 5619], [5559, 5850], [5831, 6045], [5970, 6045], [6045, 6310]]
    seq_range = seq_range +   [[6225, 8795], [8379, 8469], [8379, 8653], [8797, 9417], [9086, 9719]]
    seq_label = [       'gag',        'pol',        'vif',        'vpr', 'tat exon 1', 'rev exon 1',        'vpu']
    seq_label = seq_label +  [        'env', 'tat exon 2', 'rev exon 2',        'nef',     "3' LTR"]

    landmark_start = []
    landmark_end   = []
    landmark_label = []

    idx = int(df_index.iloc[0].HXB2)
    for i in range(len(df_index)):
        try:
            idx = int(df_index.iloc[i].HXB2)
        except:
            pass
        for j in range(len(seq_range)):
            if idx==seq_range[j][0] or (i==0 and idx>seq_range[j][0] and idx<seq_range[j][1]):
                landmark_start.append(i)
                landmark_end.append(i)
                landmark_label.append(seq_label[j])
            if idx==seq_range[j][1] or (i==len(df_index)-1 and idx>seq_range[j][0] and int(df_index.iloc[0].HXB2)<=seq_range[j][0]
                                        and idx<seq_range[j][1]):
                landmark_end[landmark_label.index(seq_label[j])] = i

    # Epitope labels

    # seq_range = epitope_range.copy()
    # seq_label = epitope_label.copy()

    # epitope_start = []
    # epitope_end   = []
    # epitope_label = []
    # epitope_sites = []
    # site2epitope  = {}

    # idx = int(df_index.iloc[0].HXB2)
    # for i in range(len(df_index)):
    #     try:
    #         idx = int(df_index.iloc[i].HXB2)
    #     except:
    #         pass
    #     for j in range(len(seq_range)):
    #         if idx==seq_range[j][0] or (i==0 and idx>seq_range[j][0] and idx<seq_range[j][1]):
    #             epitope_start.append(i)
    #             epitope_end.append(i)
    #             epitope_label.append(seq_label[j])
    #             if seq_label[j]=='DG9':
    #                 epitope_start[-1] -= 9 # account for DEP insertion
    #         if idx==seq_range[j][1] or (i==len(df_index)-1 and idx>seq_range[j][0] and int(df_index.iloc[0].HXB2)<=seq_range[j][0]
    #                                     and idx<seq_range[j][1]):
    #             epitope_end[epitope_label.index(seq_label[j])] = i
    #             iix = epitope_label.index(seq_label[j])
    #             for k in range(epitope_start[iix],epitope_end[iix]):
    #                 epitope_sites.append(k)
    #                 if k in site2epitope:
    #                     print('Unexpected! Overlapping epitopes.')
    #                 else:
    #                     site2epitope[k] = seq_label[j]
    #             print(''.join(list(df_index.iloc[epitope_start[iix]:epitope_end[iix]].TF)))


    # Dot plot for significant selection coefficients

    c_dot  = { True: C_BEN, False: C_DEL }
    s_mult = 40*SMALLSIZEDOT

    start_idx = 0
    end_idx   = len(df_index)

    def_buff = 100
    def site2angle(site, deg=False, buffer=def_buff):
        if deg: return    -360.*(site+(buffer/2))/(end_idx-start_idx+buffer)
        else:   return -2*np.pi*(site+(buffer/2))/(end_idx-start_idx+buffer)

    print(np.max(sig_s), np.min(sig_s))

    # s_max    = 0.125
    # s_min    = -0.05
    # s_norm   = s_max - s_min
    # s_levels = [-0.05, -0.025, 0, 0.025, 0.05, 0.075, 0.10, 0.125]

    # Compress s range for plotting clarity
    s_max    = 0.075
    s_min    = -0.05
    s_norm   = s_max - s_min
    s_levels = [-0.05, -0.025, 0, 0.025, 0.05, 0.075]
    sig_s    = np.array(sig_s)
    sig_s[sig_s>s_max] = s_max

    r_max = 0.90
    r_min = 0.45
    s2r = lambda s: ((s-s_min)/s_norm) * (r_max-r_min) + r_min
    level_rad = [s2r(s) for s in s_levels]

    dot_colors = [c_dot[s>0]          for s in sig_s]
    dot_sizes  = [s_mult * np.fabs(s) for s in sig_s]
    dot_rads   = [s2r(s) for s in sig_s]
    dot_x = [dot_rads[i] * np.cos(site2angle(sig_site_real[i])) for i in range(len(sig_s))]
    dot_y = [dot_rads[i] * np.sin(site2angle(sig_site_real[i])) for i in range(len(sig_s))]

    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
    s_text_levels = [-0.05, 0, 0.05]
    text_level_rad = [s2r(s) for s in s_text_levels]
    for i in range(len(s_text_levels)):
        ax.text(text_level_rad[i], 0, int(100*s_text_levels[i]), **txtprops)

    # Lines for polymorphic sites

    # line_r = [0.67, 0.71]
    # line_x = [[line_r[0] * np.cos(site2angle(i)), line_r[1] * np.cos(site2angle(i))] for i in np.unique(df.alignment_index)]
    # line_y = [[line_r[0] * np.sin(site2angle(i)), line_r[1] * np.sin(site2angle(i))] for i in np.unique(df.alignment_index)]
    # line_c = [LCOLOR                                                                 for i in np.unique(df.alignment_index)]

    # Arcs for  sequence landmarks

    arc_r            = [0.96, 0.99, 1.02]
    arc_dr           = 0.01
    track            = 0
    landmark_patches = []
    landmark_text_d  = { 'pol'    : -40,
                         'vif'    :  10,
                         "3' LTR" :  40 }

    for i in range(len(landmark_label)):
        if i>0 and landmark_start[i]<landmark_end[i-1]:
            track = (track + 1)%len(arc_r)
        
        arc_props = dict(center=[0,0], r=arc_r[track], width=arc_dr, lw=AXWIDTH/2, ec=BKCOLOR, fc='none',
                         theta1=site2angle(landmark_end[i], deg=True), theta2=site2angle(landmark_start[i], deg=True))
        landmark_patches.append(matplotlib.patches.Wedge(**arc_props))

        # label with line
        if landmark_label[i] in ['tat exon 1', 'tat exon 2', 'rev exon 1', 'rev exon 2']:
            txtprops  = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
            plotprops = dict(lw=AXWIDTH/2, ls='-', clip_on=False)

            label_x = [arc_r[track]*np.cos(site2angle(landmark_start[i]))]
            label_y = [arc_r[track]*np.sin(site2angle(landmark_start[i]))]
            ddr     = 0.11
            if landmark_label[i] in label2ddr:
                ddr = label2ddr[landmark_label[i]]

            ddx1 =  ddr * np.cos(site2angle(landmark_start[i]))
            ddy  =  ddr * np.sin(site2angle(landmark_start[i]))
            ddx2 = ddx1 + np.sign(ddx1)*0.03

            label_x = label_x + [label_x[0] + ddx1, label_x[0] + ddx2]
            label_y = label_y + [label_y[0] +  ddy, label_y[0] +  ddy]
            if label_x[0]<0:
                txtprops['ha'] = 'right'
            else:
                txtprops['ha'] = 'left'

            ax.text(label_x[-1], label_y[-1], landmark_label[i], **txtprops)
            mp.line(ax=ax, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)

        # plot normally
        else:
            delta_site = 50
            if landmark_label[i] in ['vif']:
                delta_site += landmark_text_d[landmark_label[i]] + 25
            elif landmark_label[i]=='vpu' and '703010505' in s_file:
                delta_site = 25
            elif landmark_label[i]=='env' and '703010505' in s_file:
                delta_site = 75
            elif landmark_label[i] in ["3' LTR"] or (landmark_label[i]=='pol' and '700010077' in s_file):
                delta_site += landmark_text_d[landmark_label[i]]
            txt_dr   = 0.04
            txt_ang  = site2angle(landmark_start[i]+delta_site)
            txt_rot  = site2angle(landmark_start[i]+delta_site, deg=True)-90
            txt_x    = (arc_r[track] + (arc_dr/2) + txt_dr)*np.cos(txt_ang)
            txt_y    = (arc_r[track] + (arc_dr/2) + txt_dr)*np.sin(txt_ang)
            txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY,
                            size=SIZELABEL, rotation=txt_rot)
            ax.text(txt_x, txt_y, landmark_label[i], **txtprops)

    # Arcs for TGCA selection tracks
    
    for i in range(len(level_rad)):
        arc_props = dict(center=[0,0], r=level_rad[i], width=0, lw=AXWIDTH/2, ec=C_NEU_LT, fc='none', theta1=0, theta2=360)
        landmark_patches.append(matplotlib.patches.Wedge(**arc_props))

    # Outer arcs for additional sequence annotations
        
    arc_range = [[6615, 6694], [6697, 6782], [7047, 7073], [7110, 7217], [7377, 7478], [7602, 7634]]
    arc_label = [       'V1',          'V2',     'loop D',         'V3',         'V4',         'V5']

    landmark_start = []
    landmark_end   = []
    landmark_label = []

    idx = int(df_index.iloc[0].HXB2)
    for i in range(len(df_index)):
        try:
            idx = int(df_index.iloc[i].HXB2)
        except:
            pass
        for j in range(len(arc_range)):
            if idx==arc_range[j][0] or (i==0 and idx>arc_range[j][0] and idx<arc_range[j][1]):
                landmark_start.append(i)
                landmark_end.append(i)
                landmark_label.append(arc_label[j])
            if idx==arc_range[j][1] or (i==len(df_index)-1 and idx>arc_range[j][0] and int(df_index.iloc[0].HXB2)<=arc_range[j][0]
                                        and idx<arc_range[j][1]):
                landmark_end[landmark_label.index(arc_label[j])] = i

    arc_r = 1.05
    for i in range(len(landmark_label)):
        arc_props = dict(center=[0,0], r=arc_r, width=0, lw=AXWIDTH/2, ec=BKCOLOR, fc='none',
                         theta1=site2angle(landmark_end[i], deg=True), theta2=site2angle(landmark_start[i], deg=True))
        landmark_patches.append(matplotlib.patches.Wedge(**arc_props))

        # plot normally
        delta_site = 50
        if landmark_label[i]=='V1' and '703010848' in s_file:
            delta_site = 25
        elif landmark_label[i]=='V2' and '703010848' in s_file:
            delta_site = 75
        elif landmark_label[i]=='loop D' and '703010505' in s_file:
            delta_site = 15
        elif landmark_label[i]=='loop D' and '703010848' in s_file:
            delta_site = 0
        elif landmark_label[i]=='V3' and '703010848' in s_file:
            delta_site = 65
        elif landmark_label[i]=='V5' and '703010505' in s_file:
            delta_site = 15
        elif landmark_label[i]=='V5' and '703010848' in s_file:
            delta_site = 15
        txt_dr   = 0.04
        txt_ang  = site2angle(landmark_start[i]+delta_site)
        txt_rot  = site2angle(landmark_start[i]+delta_site, deg=True)-90
        txt_x    = (arc_r + txt_dr)*np.cos(txt_ang)
        txt_y    = (arc_r + txt_dr)*np.sin(txt_ang)
        txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY,
                        size=SIZELABEL, rotation=txt_rot)
        ax.text(txt_x, txt_y, landmark_label[i], **txtprops)

    # # Arcs for epitopes

    # arc_r           = 1.04
    # arc_dr          = 0.32
    # epitope_patches = []

    # for i in range(len(epitope_label)):
        
    #     arc_props = dict(center=[0,0], r=arc_r, width=arc_dr, lw=AXWIDTH/2, ec=BKCOLOR, fc='#f2f2f2', alpha=0.8,
    #                      theta1=site2angle(epitope_end[i], deg=True), theta2=site2angle(epitope_start[i], deg=True))
    #     epitope_patches.append(matplotlib.patches.Wedge(**arc_props))

    #     # label epitopes
    #     if True:
    #         mid       = (site2angle(epitope_end[i], deg=True)+site2angle(epitope_start[i], deg=True))/2.
    #         label_x   = [arc_r * np.cos(mid * 2 * np.pi / 360.)]
    #         label_y   = [arc_r * np.sin(mid * 2 * np.pi / 360.)]
    #         ddr       = 0.06
    #         if epitope_label[i] in label2ddr:
    #             ddr = label2ddr[epitope_label[i]]
    #         ddx1 =  ddr * np.cos(mid * 2 * np.pi / 360.)
    #         ddy  =  ddr * np.sin(mid * 2 * np.pi / 360.)
    #         ddx2 = ddx1 + np.sign(ddx1)*0.03
    #         txtprops  = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
    #         plotprops = dict(lw=AXWIDTH/2, ls='-', clip_on=False)

    #         label_x = label_x + [label_x[0] + ddx1, label_x[0] + ddx2]
    #         label_y = label_y + [label_y[0] +  ddy, label_y[0] +  ddy]
    #         if label_x[0]<0:
    #             txtprops['ha'] = 'right'
    #         else:
    #             txtprops['ha'] = 'left'

    #         ax.text(label_x[-1], label_y[-1], epitope_label[i], **txtprops)
    #         mp.line(ax=ax, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)

    # # Arc plot for large values of integrated covariance

    # c_pos  = LCOLOR
    # c_neg  = LCOLOR
    # c_circ = { True : c_neg, False : c_pos }

    # arc_rad   = 0.65
    # arc_mult  = SIZELINE/2
    # arc_alpha = 0.1

    # circ_color = [c_circ[ic>0]                                     for ic in inv_cov]
    # #circ_color = [c_circ[idx_1[i] in epitope_sites or idx_2[i] in epitope_sites] for i in range(len(inv_cov))]
    # circ_rad   = [[arc_rad, arc_rad]                               for ic in inv_cov]
    # circ_arc   = [dict(lw=arc_mult * np.fabs(ic), alpha=arc_alpha) for ic in inv_cov]
    # circ_x     = [(i+(def_buff/2)) for i in idx_1]
    # circ_y     = [(i+(def_buff/2)) for i in idx_2]

    # Make plot

    pprops = { 'xlim':   [-1.05, 1.05],
               'ylim':   [-1.05, 1.05],
               'xticks': [],
               'yticks': [],
               'noaxes': True }

    plotprops = dict(lw=0, marker='o', s=dot_sizes[0], zorder=9999)
    mp.plot(type='scatter', ax=ax, x=[[dot_x[0]]], y=[[dot_y[0]]], colors=[dot_colors[0]], plotprops=plotprops, **pprops)

    for i in range(1, len(dot_x)):
        plotprops = dict(lw=0, marker='o', s=dot_sizes[i], zorder=9999)
        mp.scatter(ax=ax, x=[[dot_x[i]]], y=[[dot_y[i]]], colors=[dot_colors[i]], plotprops=plotprops)

    # for i in range(len(np.unique(df_info.alignment_index))):
    #     plotprops = dict(lw=AXWIDTH, ls='-')
    #     mp.line(ax=ax, x=[line_x[i]], y=[line_y[i]], colors=[line_c[i]], plotprops=plotprops)

    # mp.plot(type='circos', ax=ax, x=circ_x, y=circ_y, **pprops)

    for patch in landmark_patches:
        ax.add_artist(patch)

    # for patch in epitope_patches:
    #     ax.add_artist(patch)

    return sig_s, sig_site_real, sig_nuc_idx


def plot_figure_ch505_structure(filename=None):
    """
    CH505 structure and views of selection coefficients, generated in Pymol and arranged here.
    """

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN # SLIDE_WIDTH
    hshrink = 0.55
    goldh   = w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    x_side = 0.33
    x_top  = 0.28

    # hshrink = 0.6
    # x_side = 0.38
    # x_top = 0.35

    ## figure boundaries
    box_t = 0.98
    box_b = 0.02
    box_l = 0.02
    box_r = 0.98

    # midpoint boundaries
    l_end = 0.52
    r_start = 0.52

    ## side view: choose width, scale height to be square
    # x_side = 0.42
    y_side = x_side * w / goldh
    l_side = (0.5 - box_l - x_side)/2 + box_l
    box_side = dict(left=l_side, right=l_side+x_side, bottom=box_t-y_side, top=box_t)

    ## side sub views: fit two in the same space as the side view, center views along y axis
    img_aspect_ratio = 679/212
    dy_mult = 1.1
    x_side_sub = (l_end - box_l)/2
    y_side_sub = x_side_sub * (w / goldh) / img_aspect_ratio
    # y_side_sub = (box_t - y_side - box_b) / (2 * dy_mult)
    dy_side_sub = (box_t - y_side - box_b - 2*y_side_sub)/2
    box_side_sub1 = dict(left=box_l,              right=box_l+x_side_sub,       bottom=box_b+1*y_side_sub+1*dy_side_sub, top=box_b+2*y_side_sub+1*dy_side_sub)
    box_side_sub2 = dict(left=box_l+x_side_sub,   right=box_l+2*x_side_sub,     bottom=box_b+1*y_side_sub+1*dy_side_sub, top=box_b+2*y_side_sub+1*dy_side_sub)
    box_side_sub3 = dict(left=box_l+x_side_sub/2, right=(box_l+3*x_side_sub/2), bottom=box_b+0*y_side_sub+0*dy_side_sub, top=box_b+1*y_side_sub+0*dy_side_sub)

    ## top view: choose width, scale height to be square
    # x_top = 0.40
    y_top = x_top * w / goldh
    l_top = (box_r - r_start - x_top)/2 + r_start
    box_top = dict(left=l_top, right=l_top+x_top, bottom=box_t-y_top, top=box_t)

    ## top sub views: in remaining space, scale square subplots such that 3 can fit across and 2 vertically (2 top row 3 bottom row)
    dxy_mult = 1.02
    y_space = (box_t - box_b - y_top) * goldh
    x_space = (box_r - r_start) * w
    x_top_sub = 0
    y_top_sub = 0
    if x_space/3 < y_space/2:
        x_top_sub = x_space / (3 * dxy_mult * w)
        y_top_sub = x_top_sub * (w / goldh)
    else:
        y_top_sub = y_space / (2 * dxy_mult * goldh)
        x_top_sub = y_top_sub * (goldh / w)
    dx_top_sub = (box_r - r_start - 3*x_top_sub)/2
    dy_top_sub = (box_t - box_b - y_top - 2*y_top_sub)/2

    box_top_sub1 = dict(left=r_start+0.5*dx_top_sub+0.5*x_top_sub, right=r_start+0.5*dx_top_sub+1.5*x_top_sub, bottom=box_b+1*dy_top_sub+1*y_top_sub, top=box_b+1*dy_top_sub+2*y_top_sub)
    box_top_sub2 = dict(left=r_start+1.5*dx_top_sub+1.5*x_top_sub, right=r_start+1.5*dx_top_sub+2.5*x_top_sub, bottom=box_b+1*dy_top_sub+1*y_top_sub, top=box_b+1*dy_top_sub+2*y_top_sub)
    box_top_sub3 = dict(left=r_start+0*dx_top_sub+0*x_top_sub,     right=r_start+0*dx_top_sub+1*x_top_sub,     bottom=box_b+0*dy_top_sub+0*y_top_sub, top=box_b+0*dy_top_sub+1*y_top_sub)
    box_top_sub4 = dict(left=r_start+1*dx_top_sub+1*x_top_sub,     right=r_start+1*dx_top_sub+2*x_top_sub,     bottom=box_b+0*dy_top_sub+0*y_top_sub, top=box_b+0*dy_top_sub+1*y_top_sub)
    box_top_sub5 = dict(left=r_start+2*dx_top_sub+2*x_top_sub,     right=r_start+2*dx_top_sub+3*x_top_sub,     bottom=box_b+0*dy_top_sub+0*y_top_sub, top=box_b+0*dy_top_sub+1*y_top_sub)

    ## create gridspecs and axes
    gs_side      = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_side)
    gs_side_sub1 = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_side_sub1)
    gs_side_sub2 = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_side_sub2)
    gs_side_sub3 = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_side_sub3)
    gs_top       = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_top)
    gs_top_sub1  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_top_sub1)
    gs_top_sub2  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_top_sub2)
    gs_top_sub3  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_top_sub3)
    gs_top_sub4  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_top_sub4)
    gs_top_sub5  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_top_sub5)

    ax_side      = plt.subplot(gs_side[0, 0])
    ax_side_sub1 = plt.subplot(gs_side_sub1[0, 0])
    ax_side_sub2 = plt.subplot(gs_side_sub2[0, 0])
    ax_side_sub3 = plt.subplot(gs_side_sub3[0, 0])
    ax_top       = plt.subplot(gs_top[0, 0])
    ax_top_sub1  = plt.subplot(gs_top_sub1[0, 0])
    ax_top_sub2  = plt.subplot(gs_top_sub2[0, 0])
    ax_top_sub3  = plt.subplot(gs_top_sub3[0, 0])
    ax_top_sub4  = plt.subplot(gs_top_sub4[0, 0])
    ax_top_sub5  = plt.subplot(gs_top_sub5[0, 0])

    # PLOT IMAGES

    struct_dir = 'data/structure/CH505'

    def make_square(img, trim_top=13, trim_bottom=5, set_min_w=None, pad=False):
        img = img[trim_top:-trim_bottom, :, :]
        h   = img.shape[0]
        mid = img.shape[1]//2
        low_lim = mid - h//2
        high_lim = mid + h//2
        if set_min_w is not None:
            if low_lim>set_min_w[0]:
                low_lim = set_min_w[0]
            if high_lim<set_min_w[1]:
                high_lim = set_min_w[1]
        img = img[:, low_lim:high_lim, :]
        # pad vertically so that the image is square
        if pad:
            h   = img.shape[0]
            w   = img.shape[1]
            pad = (w-h)//2
            img = np.pad(img, ((0,0), (pad,pad), (0,0)), mode='constant')
        return img

    def make_landscape(img, center_y=0.25, trim_top=5, trim_bottom=5, aspect_ratio=(y_side_sub*goldh)/(x_side_sub*w), set_min_w=[320, 1140]):
        img = img[trim_top:-trim_bottom, :, :]
        if set_min_w is not None:
            img = img[:, set_min_w[0]:set_min_w[1], :]
        mid = int(img.shape[0]*center_y)
        dh  = int(img.shape[1]*aspect_ratio)
        low_lim = np.max([0, mid-dh//2])
        high_lim = np.min([img.shape[0], mid+dh//2])
        return img[low_lim:high_lim, :, :]
    
    def trim_to_size(img, trim_top=0, trim_bottom=0, trim_left=0, trim_right=0):
        if trim_top>0 or trim_bottom>0:
            img = img[trim_top:-trim_bottom, :, :]
        if trim_left>0 or trim_right>0:
            img = img[:, trim_left:-trim_right, :]
        return img
    
    tprops = dict(family=FONTFAMILY, size=SIZELABEL, color=BKCOLOR, ha='center', va='center')
    
    side_trim_top = 85
    img = mpimg.imread('%s/side-view/entire-view.png' % struct_dir)
    img = make_square(img, trim_top=side_trim_top, set_min_w=[358, 1041], pad=True)
    ax_side.imshow(img)

    img = mpimg.imread('%s/side-view/only-CH103_clipped.png' % struct_dir)
    img = trim_to_size(img, 0, 0, 0, 0)
    # img = make_landscape(img, trim_top=side_trim_top)
    ax_side_sub1.imshow(img)
    ax_side_sub1.text((box_side_sub1['left']+box_side_sub1['right'])/2, box_side_sub1['top']+0.01, 'CH103 binding site', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/side-view/only-CH235_clipped.png' % struct_dir)
    img = trim_to_size(img, 1, 1, 1, 1)
    # img = make_landscape(img, trim_top=side_trim_top)
    ax_side_sub2.imshow(img)
    ax_side_sub2.text((box_side_sub2['left']+box_side_sub2['right'])/2, box_side_sub2['top']+0.01, 'CH235 binding site', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/side-view/only-autologous_clipped.png' % struct_dir)
    img = trim_to_size(img, 2, 1, 0, 0)
    # img = make_landscape(img, trim_top=side_trim_top)
    ax_side_sub3.imshow(img)
    ax_side_sub3.text((box_side_sub3['left']+box_side_sub3['right'])/2, box_side_sub3['top']+0.01, 'Strain-specific Ab escape', transform=fig.transFigure, **tprops)

    top_trim_bottom = 50
    img = mpimg.imread('%s/top-view/entire-view.png' % struct_dir)
    img = make_square(img, trim_bottom=top_trim_bottom)
    ax_top.imshow(img)

    img = mpimg.imread('%s/top-view/only-V1V2.png' % struct_dir)
    img = make_square(img, trim_bottom=top_trim_bottom)
    ax_top_sub1.imshow(img)
    ax_top_sub1.text((box_top_sub1['left']+box_top_sub1['right'])/2, box_top_sub1['top']+0.01, 'V1-V2', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/top-view/only-V3.png' % struct_dir)
    img = make_square(img, trim_bottom=top_trim_bottom)
    ax_top_sub2.imshow(img)
    ax_top_sub2.text((box_top_sub2['left']+box_top_sub2['right'])/2, box_top_sub2['top']+0.01, 'V3', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/top-view/only-V4.png' % struct_dir)
    img = make_square(img, trim_bottom=top_trim_bottom)
    ax_top_sub3.imshow(img)
    ax_top_sub3.text((box_top_sub3['left']+box_top_sub3['right'])/2, box_top_sub3['top']+0.01, 'V4', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/top-view/only-V5.png' % struct_dir)
    img = make_square(img, trim_bottom=top_trim_bottom)
    ax_top_sub4.imshow(img)
    ax_top_sub4.text((box_top_sub4['left']+box_top_sub4['right'])/2, box_top_sub4['top']+0.01, 'V5', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/top-view/only-CD4BS.png' % struct_dir)
    img = make_square(img, trim_bottom=top_trim_bottom)
    ax_top_sub5.imshow(img)
    ax_top_sub5.text((box_top_sub5['left']+box_top_sub5['right'])/2, box_top_sub5['top']+0.01, 'CD4 binding site', transform=fig.transFigure, **tprops)

    # Add plot sublabels
    dy_sub = 0.02
    dx_sub = -0.02

    ax_side.text(box_side['left']+dy_sub, box_side['top']+dx_sub, 'a'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_side_sub1.text(box_side_sub1['left']+dy_sub, box_side_sub1['top']+dy_sub, 'b'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_side_sub2.text(box_side_sub2['left']+dy_sub, box_side_sub2['top']+dy_sub, 'c'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_side_sub3.text(box_side_sub3['left']+dy_sub, box_side_sub3['top']+dy_sub-0.01, 'd'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ax_top.text(box_top['left']+dy_sub, box_top['top']+dx_sub, 'e'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_top_sub1.text(box_top_sub1['left']+dy_sub, box_top_sub1['top']+dx_sub, 'f'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_top_sub2.text(box_top_sub2['left']+dy_sub, box_top_sub2['top']+dx_sub, 'g'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_top_sub3.text(box_top_sub3['left']+dy_sub, box_top_sub3['top']+dx_sub, 'h'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_top_sub4.text(box_top_sub4['left']+dy_sub, box_top_sub4['top']+dx_sub, 'i'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_top_sub5.text(box_top_sub5['left']+dy_sub, box_top_sub5['top']+dx_sub, 'j'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    # REMOVE TICKS AND SPINES

    for ax in [ax_side, ax_side_sub1, ax_side_sub2, ax_side_sub3, ax_top, ax_top_sub1, ax_top_sub2, ax_top_sub3, ax_top_sub4, ax_top_sub5]:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

    # Add legend
    
    print(ax_side.get_xlim(), ax_side.get_ylim())

    # ## color bar
    ddx      = 0 #0.525
    ddy      = 37
    label_x  = 810
    label_y  = 490
    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_side.text(label_x, label_y+ddy, 'Inferred selection\ncoefficient (%)', **txtprops)
    
    aspect_ratio = (ax_side.get_xlim()[1]-ax_side.get_xlim()[0]) / (ax_side.get_ylim()[0]-ax_side.get_ylim()[1])
    dbox      = 20
    rec_props = dict(height=dbox, width=dbox, ec=None, lw=AXWIDTH/2, clip_on=False)
    for i in range(-5, 5+1, 1):
        c = BKCOLOR
        t = i/5
        if t>0:
            c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
        else:
            c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
        rec = matplotlib.patches.Rectangle(xy=((i-0.5)*dbox + label_x, label_y + 2*ddy), fc=c, **rec_props)
        ax_side.add_artist(rec)

    txtprops = dict(ha='center', va='top', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_side.text(-5.25*dbox + label_x, label_y + 3*ddy, -8, **txtprops)
    ax_side.text(             label_x, label_y + 3*ddy,  0, **txtprops)
    ax_side.text( 5.00*dbox + label_x, label_y + 3*ddy,  8, **txtprops)

#     txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
#     for i in range(len(NUC)-1):
#         ax_smpl.text(-0.85, 3-i+0.5, NUC[i+1], clip_on=False, **txtprops)

#     txtprops['ha'] = 'left'
#     ax_smpl.text(1.3, -3.5, 'TF nucleotide', clip_on=False, **txtprops)
#     ax_smpl.text(1.3, -3.5 + 4*coef_legend_dy, 'Not observed', clip_on=False, **txtprops) # paper
# #    ax_smpl.text(1.3, -3.5 + 6*coef_legend_dy, 'Not observed', clip_on=False, **txtprops) # slides

#     txtprops['ha'] = 'center'
#     txtprops['va'] = 'top'
#     for i in range((epitope_end[eidx]-epitope_start[eidx]+1)//3):
#         ax_smpl.text(1.5+i*sub_dx, -0.5, 160+i, clip_on=False, **txtprops)

#     ax_smpl.text(sub_box*sub_dx - 11, -4.5, -5, clip_on=False, **txtprops)
#     ax_smpl.text(sub_box*sub_dx -  6, -4.5,  0, clip_on=False, **txtprops)
#     ax_smpl.text(sub_box*sub_dx -  1, -4.5,  5, clip_on=False, **txtprops)
#     ax_smpl.text(sub_box*sub_dx -  5.5, -6.0, 'Inferred selection\ncoefficient, $\hat{s}$ (%)', clip_on=False, **txtprops)

    # SAVE FIGURE

    if filename==None:
        filename = 'fig-ch505-structure%s' % EXT

    plt.savefig('%s%s' % (FIG_DIR, filename), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('CH505-structure done.')


def plot_figure_ch848_structure(filename=None):
    """
    CH848 structure and views of selection coefficients, generated in Pymol and arranged here.
    """

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN # SLIDE_WIDTH
    hshrink = 0.55
    goldh   = w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    x_side = 0.33
    x_top  = 0.28

    # hshrink = 0.6
    # x_side = 0.38
    # x_top = 0.35

    ## figure boundaries
    box_t = 0.98
    box_b = 0.02
    box_l = 0.02
    box_r = 0.98

    # midpoint boundaries
    l_end = 0.52
    r_start = 0.52

    ## side view: choose width, scale height to be square
    # x_side = 0.42
    y_side = x_side * w / goldh
    l_side = (0.5 - box_l - x_side)/2 + box_l
    box_side = dict(left=l_side, right=l_side+x_side, bottom=box_t-y_side+0.02, top=box_t+0.02)

    ## side sub views: two by two square subplots
    img_aspect_ratio = 1
    x_side_sub = (l_end - box_l)/4.5
    y_side_sub = x_side_sub * (w / goldh) / img_aspect_ratio
    # y_side_sub = (box_t - y_side - box_b) / (2 * dy_mult)
    dy_side_sub = (box_t - y_side - box_b - 2*y_side_sub)/2
    box_l = (l_side + l_side + x_side)/2 - 1.1*x_side_sub
    box_side_sub1 = dict(left=box_l,                right=box_l+x_side_sub,   bottom=box_b+1*y_side_sub+1*dy_side_sub+0.04, top=box_b+2*y_side_sub+1*dy_side_sub+0.04)
    box_side_sub2 = dict(left=box_l+1.2*x_side_sub, right=box_l+2.2*x_side_sub, bottom=box_b+1*y_side_sub+1*dy_side_sub+0.04, top=box_b+2*y_side_sub+1*dy_side_sub+0.04)
    box_side_sub3 = dict(left=box_l,                right=box_l+x_side_sub,   bottom=box_b+0*y_side_sub+0*dy_side_sub, top=box_b+1*y_side_sub+0*dy_side_sub)
    box_side_sub4 = dict(left=box_l+1.2*x_side_sub, right=box_l+2.2*x_side_sub, bottom=box_b+0*y_side_sub+0*dy_side_sub, top=box_b+1*y_side_sub+0*dy_side_sub)

    ## top view: choose width, scale height to be square
    # x_top = 0.40
    y_top = x_top * w / goldh
    l_top = (box_r - r_start - x_top)/2 + r_start
    box_top = dict(left=l_top, right=l_top+x_top, bottom=box_t-y_top, top=box_t)

    ## top sub views: in remaining space, scale square subplots such that 3 can fit across and 2 vertically (2 top row 3 bottom row)
    dxy_mult = 1.02
    y_space = (box_t - box_b - y_top) * goldh
    x_space = (box_r - r_start) * w
    x_top_sub = 0
    y_top_sub = 0
    if x_space/3 < y_space/2:
        x_top_sub = x_space / (3 * dxy_mult * w)
        y_top_sub = x_top_sub * (w / goldh)
    else:
        y_top_sub = y_space / (2 * dxy_mult * goldh)
        x_top_sub = y_top_sub * (goldh / w)
    dx_top_sub = (box_r - r_start - 3*x_top_sub)/2
    dy_top_sub = (box_t - box_b - y_top - 2*y_top_sub)/2

    box_top_sub1 = dict(left=r_start+0.5*dx_top_sub+0.5*x_top_sub, right=r_start+0.5*dx_top_sub+1.5*x_top_sub, bottom=box_b+1*dy_top_sub+1*y_top_sub, top=box_b+1*dy_top_sub+2*y_top_sub)
    box_top_sub2 = dict(left=r_start+1.5*dx_top_sub+1.5*x_top_sub, right=r_start+1.5*dx_top_sub+2.5*x_top_sub, bottom=box_b+1*dy_top_sub+1*y_top_sub, top=box_b+1*dy_top_sub+2*y_top_sub)
    box_top_sub3 = dict(left=r_start+0*dx_top_sub+0*x_top_sub,     right=r_start+0*dx_top_sub+1*x_top_sub,     bottom=box_b+0*dy_top_sub+0*y_top_sub, top=box_b+0*dy_top_sub+1*y_top_sub)
    box_top_sub4 = dict(left=r_start+1*dx_top_sub+1*x_top_sub,     right=r_start+1*dx_top_sub+2*x_top_sub,     bottom=box_b+0*dy_top_sub+0*y_top_sub, top=box_b+0*dy_top_sub+1*y_top_sub)
    box_top_sub5 = dict(left=r_start+2*dx_top_sub+2*x_top_sub,     right=r_start+2*dx_top_sub+3*x_top_sub,     bottom=box_b+0*dy_top_sub+0*y_top_sub, top=box_b+0*dy_top_sub+1*y_top_sub)

    ## create gridspecs and axes
    gs_side      = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_side)
    gs_side_sub1 = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_side_sub1)
    gs_side_sub2 = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_side_sub2)
    gs_side_sub3 = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_side_sub3)
    gs_side_sub4 = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_side_sub4)
    gs_top       = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_top)
    gs_top_sub1  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_top_sub1)
    gs_top_sub2  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_top_sub2)
    gs_top_sub3  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_top_sub3)
    gs_top_sub4  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_top_sub4)
    gs_top_sub5  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_top_sub5)

    ax_side      = plt.subplot(gs_side[0, 0])
    ax_side_sub1 = plt.subplot(gs_side_sub1[0, 0])
    ax_side_sub2 = plt.subplot(gs_side_sub2[0, 0])
    ax_side_sub3 = plt.subplot(gs_side_sub3[0, 0])
    ax_side_sub4 = plt.subplot(gs_side_sub4[0, 0])
    ax_top       = plt.subplot(gs_top[0, 0])
    ax_top_sub1  = plt.subplot(gs_top_sub1[0, 0])
    ax_top_sub2  = plt.subplot(gs_top_sub2[0, 0])
    ax_top_sub3  = plt.subplot(gs_top_sub3[0, 0])
    ax_top_sub4  = plt.subplot(gs_top_sub4[0, 0])
    ax_top_sub5  = plt.subplot(gs_top_sub5[0, 0])

    # PLOT IMAGES

    struct_dir = 'data/structure/CH848'

    def make_square(img, trim_top=13, trim_bottom=5, set_min_w=None, pad=False):
        img = img[trim_top:-trim_bottom, :, :]
        h   = img.shape[0]
        mid = img.shape[1]//2
        low_lim = mid - h//2
        high_lim = mid + h//2
        if set_min_w is not None:
            if low_lim>set_min_w[0]:
                low_lim = set_min_w[0]
            if high_lim<set_min_w[1]:
                high_lim = set_min_w[1]
        img = img[:, low_lim:high_lim, :]
        # pad vertically so that the image is square
        if pad:
            h   = img.shape[0]
            w   = img.shape[1]
            pad = (w-h)//2
            img = np.pad(img, ((0,0), (pad,pad), (0,0)), mode='constant')
        return img

    def make_landscape(img, center_y=0.25, trim_top=5, trim_bottom=5, aspect_ratio=(y_side_sub*goldh)/(x_side_sub*w), set_min_w=[320, 1140]):
        img = img[trim_top:-trim_bottom, :, :]
        if set_min_w is not None:
            img = img[:, set_min_w[0]:set_min_w[1], :]
        mid = int(img.shape[0]*center_y)
        dh  = int(img.shape[1]*aspect_ratio)
        low_lim = np.max([0, mid-dh//2])
        high_lim = np.min([img.shape[0], mid+dh//2])
        return img[low_lim:high_lim, :, :]
    
    def trim_to_size(img, trim_top=0, trim_bottom=0, trim_left=0, trim_right=0):
        if trim_top>0 or trim_bottom>0:
            img = img[trim_top:-trim_bottom, :, :]
        if trim_left>0 or trim_right>0:
            img = img[:, trim_left:-trim_right, :]
        return img
    
    tprops = dict(family=FONTFAMILY, size=SIZELABEL, color=BKCOLOR, ha='center', va='center')
    
    side_trim_top = 80
    side_trim_bottom = 1
    img = mpimg.imread('%s/side-view/entire-view.png' % struct_dir)
    # print('side-view/entire-view.png', img.shape)
    img = make_square(img, trim_top=side_trim_top, trim_bottom=side_trim_bottom, set_min_w=[320, 1000], pad=True)
    ax_side.imshow(img)

    trim_left = 650
    trim_right = 230
    trim_top = 65
    trim_bottom = 374

    img = mpimg.imread('%s/side-view/only-DH270BS.png' % struct_dir)
    # print('side-view/only-DH270BS.png', img.shape)
    img = trim_to_size(img, trim_top, trim_bottom, trim_left, trim_right)
    # img = make_landscape(img, trim_top=side_trim_top)
    # img = make_square(img, trim_top=side_trim_top, trim_bottom=side_trim_bottom, set_min_w=[700, 800], pad=True)
    ax_side_sub1.imshow(img)
    ax_side_sub1.text((box_side_sub1['left']+box_side_sub1['right'])/2, box_side_sub1['top']-0.02, 'DH270', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/side-view/only-DH272BS.png' % struct_dir)
    # print('side-view/only-DH272BS.png', img.shape)
    img = trim_to_size(img, trim_top, trim_bottom, trim_left, trim_right)
    # img = make_landscape(img, trim_top=side_trim_top)
    # img = make_square(img, trim_top=side_trim_top, trim_bottom=side_trim_bottom, set_min_w=[320, 1000], pad=True)
    ax_side_sub2.imshow(img)
    ax_side_sub2.text((box_side_sub2['left']+box_side_sub2['right'])/2, box_side_sub2['top']-0.02, 'DH272', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/side-view/only-DH475BS.png' % struct_dir)
    # print('side-view/only-DH475BS.png', img.shape)
    img = trim_to_size(img, trim_top, trim_bottom, trim_left, trim_right)
    # img = make_landscape(img, trim_top=side_trim_top)
    # img = make_square(img, trim_top=side_trim_top, trim_bottom=side_trim_bottom, set_min_w=[320, 1000], pad=True)
    ax_side_sub3.imshow(img)
    ax_side_sub3.text((box_side_sub3['left']+box_side_sub3['right'])/2, box_side_sub3['top']-0.02, 'DH475', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/side-view/only-autologous.png' % struct_dir)
    # print('side-view/only-autologous.png', img.shape)
    img = trim_to_size(img, trim_top+3+10, trim_bottom+3-10, trim_left+3-10, trim_right+3+10)
    # img = make_landscape(img, trim_top=side_trim_top)
    # img = make_square(img, trim_top=side_trim_top, trim_bottom=side_trim_bottom, set_min_w=[320, 1000], pad=True)
    ax_side_sub4.imshow(img)
    ax_side_sub4.text((box_side_sub4['left']+box_side_sub4['right'])/2, box_side_sub4['top']-0.02, 'Strain-specific Ab', transform=fig.transFigure, **tprops)

    top_trim_bottom = 50
    img = mpimg.imread('%s/top-view/entire-view.png' % struct_dir)
    img = make_square(img, trim_bottom=top_trim_bottom)
    ax_top.imshow(img)

    img = mpimg.imread('%s/top-view/only-V1-and-V2.png' % struct_dir)
    img = make_square(img, trim_bottom=top_trim_bottom)
    ax_top_sub1.imshow(img)
    ax_top_sub1.text((box_top_sub1['left']+box_top_sub1['right'])/2, box_top_sub1['top']+0.01, 'V1-V2', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/top-view/only-V3.png' % struct_dir)
    img = make_square(img, trim_bottom=top_trim_bottom)
    ax_top_sub2.imshow(img)
    ax_top_sub2.text((box_top_sub2['left']+box_top_sub2['right'])/2, box_top_sub2['top']+0.01, 'V3', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/top-view/only-V4.png' % struct_dir)
    img = make_square(img, trim_bottom=top_trim_bottom)
    ax_top_sub3.imshow(img)
    ax_top_sub3.text((box_top_sub3['left']+box_top_sub3['right'])/2, box_top_sub3['top']+0.01, 'V4', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/top-view/only-V5.png' % struct_dir)
    img = make_square(img, trim_bottom=top_trim_bottom)
    ax_top_sub4.imshow(img)
    ax_top_sub4.text((box_top_sub4['left']+box_top_sub4['right'])/2, box_top_sub4['top']+0.01, 'V5', transform=fig.transFigure, **tprops)

    img = mpimg.imread('%s/top-view/only-CD4BS.png' % struct_dir)
    img = make_square(img, trim_bottom=top_trim_bottom)
    ax_top_sub5.imshow(img)
    ax_top_sub5.text((box_top_sub5['left']+box_top_sub5['right'])/2, box_top_sub5['top']+0.01, 'CD4 binding site', transform=fig.transFigure, **tprops)

    # Add plot sublabels
    dy_sub = 0.02
    dx_sub = -0.02

    ax_side.text(box_side['left']+dy_sub, box_side['top']+dx_sub-0.02, 'a'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_side_sub1.text(box_side_sub1['left']-0.01, box_side_sub1['top']-0.02, 'b'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_side_sub2.text(box_side_sub2['left']-0.01, box_side_sub2['top']-0.02, 'c'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_side_sub3.text(box_side_sub3['left']-0.01, box_side_sub3['top']-0.02, 'd'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_side_sub4.text(box_side_sub4['left']-0.01, box_side_sub4['top']-0.02, 'e'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ax_top.text(box_top['left']+dy_sub, box_top['top']+dx_sub, 'f'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_top_sub1.text(box_top_sub1['left']+dy_sub, box_top_sub1['top']+dx_sub, 'g'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_top_sub2.text(box_top_sub2['left']+dy_sub, box_top_sub2['top']+dx_sub, 'h'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_top_sub3.text(box_top_sub3['left']+dy_sub, box_top_sub3['top']+dx_sub, 'i'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_top_sub4.text(box_top_sub4['left']+dy_sub, box_top_sub4['top']+dx_sub, 'j'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_top_sub5.text(box_top_sub5['left']+dy_sub, box_top_sub5['top']+dx_sub, 'k'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    # REMOVE TICKS AND SPINES

    for ax in [ax_side, ax_side_sub1, ax_side_sub2, ax_side_sub3, ax_side_sub4, ax_top, ax_top_sub1, ax_top_sub2, ax_top_sub3, ax_top_sub4, ax_top_sub5]:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

    # Add legend
    
    print(ax_side.get_xlim(), ax_side.get_ylim())

    # ## color bar
    ddx      = 0 #0.525
    ddy      = 37*1.1
    label_x  = 810
    label_y  = 490
    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_side.text(label_x, label_y+ddy, 'Inferred selection\ncoefficient (%)', **txtprops)
    
    aspect_ratio = (ax_side.get_xlim()[1]-ax_side.get_xlim()[0]) / (ax_side.get_ylim()[0]-ax_side.get_ylim()[1])
    dbox      = 20*1.1
    rec_props = dict(height=dbox, width=dbox, ec=None, lw=AXWIDTH/2, clip_on=False)
    for i in range(-5, 5+1, 1):
        c = BKCOLOR
        t = i/5
        if t>0:
            c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
        else:
            c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
        rec = matplotlib.patches.Rectangle(xy=((i-0.5)*dbox + label_x, label_y + 2*ddy), fc=c, **rec_props)
        ax_side.add_artist(rec)

    txtprops = dict(ha='center', va='top', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_side.text(-5.25*dbox + label_x, label_y + 3*ddy, -8, **txtprops)
    ax_side.text(             label_x, label_y + 3*ddy,  0, **txtprops)
    ax_side.text( 5.00*dbox + label_x, label_y + 3*ddy,  8, **txtprops)

#     txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
#     for i in range(len(NUC)-1):
#         ax_smpl.text(-0.85, 3-i+0.5, NUC[i+1], clip_on=False, **txtprops)

#     txtprops['ha'] = 'left'
#     ax_smpl.text(1.3, -3.5, 'TF nucleotide', clip_on=False, **txtprops)
#     ax_smpl.text(1.3, -3.5 + 4*coef_legend_dy, 'Not observed', clip_on=False, **txtprops) # paper
# #    ax_smpl.text(1.3, -3.5 + 6*coef_legend_dy, 'Not observed', clip_on=False, **txtprops) # slides

#     txtprops['ha'] = 'center'
#     txtprops['va'] = 'top'
#     for i in range((epitope_end[eidx]-epitope_start[eidx]+1)//3):
#         ax_smpl.text(1.5+i*sub_dx, -0.5, 160+i, clip_on=False, **txtprops)

#     ax_smpl.text(sub_box*sub_dx - 11, -4.5, -5, clip_on=False, **txtprops)
#     ax_smpl.text(sub_box*sub_dx -  6, -4.5,  0, clip_on=False, **txtprops)
#     ax_smpl.text(sub_box*sub_dx -  1, -4.5,  5, clip_on=False, **txtprops)
#     ax_smpl.text(sub_box*sub_dx -  5.5, -6.0, 'Inferred selection\ncoefficient, $\hat{s}$ (%)', clip_on=False, **txtprops)

    # SAVE FIGURE

    if filename==None:
        filename = 'fig-ch848-structure%s' % EXT

    plt.savefig('%s%s' % (FIG_DIR, filename), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('CH848-structure done.')


def plot_trajectory_selection(s_files, traj_files, tags, filename='fig-s-vs-time'):
    """ Plot inferred selection coefficients and times they were first observed, together with selected trajectories. """
    
    # PLOT FIGURE
    
    ## set up figure grid
    
    w     = SINGLE_COLUMN #SLIDE_WIDTH
    goldh = w / 1.35
    fig   = plt.figure(figsize=(w, goldh))

    left_bd   = 0.12
    right_bd  = 0.95
    x_space   = 0.15
    box_x     = (right_bd - left_bd - 1*x_space) * (3/4)
    top_bd    = 0.95
    bottom_bd = 0.15
    y_space   = 0.20
    box_y     = (top_bd - bottom_bd - 1*y_space)/2
    
    box_s    = [dict(left=left_bd + box_x + x_space, right=right_bd, top=top_bd - k*box_y - k*y_space, bottom=top_bd - (k+1)*box_y - k*y_space) for k in range(2)]
    box_traj = [dict(left=left_bd, right=left_bd + box_x, top=top_bd - k*box_y - k*y_space, bottom=top_bd - (k+1)*box_y - k*y_space) for k in range(2)]
    
    ## a -- selection coefficients and frequency trajectories
    
    gs_s = [gridspec.GridSpec(1, 1, **box_s[k]) for k in range(len(tags))]
    ax_s = [plt.subplot(gs_s[k][0, 0]) for k in range(len(tags))]

    gs_traj = [gridspec.GridSpec(1, 1, **box_traj[k]) for k in range(len(tags))]
    ax_traj = [plt.subplot(gs_traj[k][0, 0]) for k in range(len(tags))]
    
    s_ylim        = [-0.03, 0.10]
    s_yticks      = [-0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.10]
    s_yticklabels = [int(i) for i in 100*np.array(s_yticks)]

    s_xlim = [0, 6]

    t_ylim        = [0, 1.05]
    t_yticks      = [0, 0.5, 1]
    t_yticklabels = [0, 50, 100]

    t_xlim   = [0, 1850]
    t_xticks = [0, 250, 500, 750, 1000, 1250, 1500, 1750]

    ab_palette  = sns.color_palette(palette='flare', n_colors=5)
    alt_palette = sns.husl_palette(3)
    
    for k in range(len(tags)):

        # load data
        df_s = pd.read_csv(s_files[k])
        df_traj = pd.read_csv(traj_files[k])

        ids     = np.unique(df_traj['id'])
        cat     = []
        c_color = []
        c_label = []
        if tags[k]=='CH505':
            cat     = ['CTL', 'CH103', 'CH235', 'autologous', 'glycan']
            c_color = [C_MPL, ab_palette[1], ab_palette[3], alt_palette[1], alt_palette[2]]
            c_label = ['CTL', 'CH103', 'CH235',       'ssAb', 'glycan']
        elif tags[k]=='CH848':
            cat     = ['DH270', 'DH272', 'DH475', 'autologous', 'glycan']
            c_color = [ab_palette[0], ab_palette[2], ab_palette[4], alt_palette[1], alt_palette[2]]
            c_label = ['DH270', 'DH272', 'DH475',       'ssAb', 'glycan']

        x_std  = 0.1
        x_vals = []
        s_vals = []
        for c in cat:
            s_vals.append(df_s[df_s['category']==c]['selection'])
            x_vals.append(np.random.normal(cat.index(c)+0.5, x_std, len(s_vals[-1])))
    
        ## a/b/c.1 -- frequency trajectory
        
        times  = []
        freqs  = []
        colors = []
        
        t_smooth = 100
        for id in ids:
            # temp_times = np.array(df_traj[df_traj['id']==id]['time'])
            # temp_freqs = np.array(df_traj[df_traj['id']==id]['frequency'])

            # x = np.linspace(np.min(temp_times), np.max(temp_times), 500)
            # y = [np.sum(temp_freqs * np.exp(-np.fabs(t - temp_times)/t_smooth) / np.sum(np.exp(-np.fabs(t - temp_times)/t_smooth))) for t in x]

            # times.append(x)
            # freqs.append(y)

            times.append(df_traj[df_traj['id']==id]['time'])
            freqs.append(df_traj[df_traj['id']==id]['frequency'])
            colors.append(c_color[cat.index(id.split('.')[0])])
        
        lineprops = {'lw': SIZELINE, 'ls': '-', 'alpha': 0.5 }
        
        pprops = { 'xlim':        [lim for lim in t_xlim],
                   'xticks':      [xt for xt in t_xticks],
                   'ylim':        [lim for lim in t_ylim],
                   'yticks':      [yt for yt in t_yticks],
                   'yticklabels': [yl for yl in t_yticklabels],
                   'yminorticks': [0.25, 0.75],
                   'xlabel':      'Time (days after infection)',
                   'ylabel':      'Frequency (%)',
                   'axoffset':    0.1,
                   'theme':       'open' }

        mp.plot(type='line', ax=ax_traj[k], x=times, y=freqs, colors=colors, plotprops=lineprops, **pprops)
        
        ## a/b/c.2 -- inferred selection coefficient

        scatterprops = dict(lw=0, s=SMALLSIZEDOT, marker='o', alpha=0.25)
        errorprops   = dict(lw=AXWIDTH, markersize=SMALLSIZEDOT, fmt='.', elinewidth=AXWIDTH, capthick=0, capsize=0, mew=AXWIDTH, alpha=1)

        pprops = { 'xlim':        [lim for lim in s_xlim],
                   'xticks':      np.array(range(len(cat)))+0.5,
                   'xticklabels': c_label,
                   'ylim':        [lim for lim in s_ylim],
                   'yticks':      [yt for yt in s_yticks],
                   'yticklabels': [yl for yl in s_yticklabels],
                #    'yminorticks': yminorticks,
                   'ylabel':      'Inferred selection\ncoefficients (%)',
                   'axoffset':    0.1,
                   'theme':       'open' }

        mp.scatter(ax=ax_s[k], x=x_vals, y=s_vals, colors=c_color, plotprops=scatterprops, **pprops)

        
        mp.plot(type='error', ax=ax_s[k], x=[[i+0.5] for i in range(len(c_color))], y=[[np.mean(sv)] for sv in s_vals], yerr=[[np.std(sv)] for sv in s_vals],
                edgecolor=[BKCOLOR for i in c_color], facecolor=[i for i in c_color], plotprops=errorprops, **pprops)
        
        for tick in ax_s[k].get_xticklabels():
            tick.set_rotation(90)

    ## legends and labels
    
    sublabels = ['a', 'b']
    ppts = ['CH505', 'CH848']
    dx = -0.08
    dy =  0.02
    for k in range(len(sublabels)):
        ax_traj[k].text(box_traj[k]['left'] + dx, box_traj[k]['top'] + dy, sublabels[k].upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
        ax_traj[k].text(box_traj[k]['left'], box_traj[k]['top'], ppts[k], transform=fig.transFigure, **DEF_LABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, filename, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)


def plot_trajectory_expanded(traj_files, tags, filename='fig-trajectory-expanded'):
    """ Plot expanded frequency trajectories over time. """
    
    # PLOT FIGURE
    
    ## set up figure grid

    w     = DOUBLE_COLUMN #SLIDE_WIDTH
    goldh = w / 1.5
    fig   = plt.figure(figsize=(w, goldh))

    top_bd    = 0.97
    bottom_bd = 0.08
    x_space   = 0.10
    left_bd   = 0.10
    right_bd  = 0.95
    x_size    = (right_bd - left_bd - x_space)/2
    
    box_traj = [dict(left=left_bd+k*(x_size+x_space), right=left_bd+k*(x_size+x_space)+x_size, top=top_bd, bottom=bottom_bd) for k in range(2)]
    
    ## a -- frequency over time by category

    categories = ['BnAb lineage resistance', 'Strain-specific Ab resistance', 'N-linked glycosylation site', 'CTL epitope']
    tag2catidx = dict(CH103=0, CH235=0, DH270=0, DH272=0, DH475=0, autologous=1, glycan=2, CTL=3)

    ab_palette   = sns.color_palette(palette='flare', n_colors=5)
    alt_palette  = sns.husl_palette(3)
    catidx2color = [ab_palette[1], alt_palette[1], alt_palette[2], C_MPL]
    
    gs_traj = [gridspec.GridSpec(len(categories), 1, **box_traj[k]) for k in range(len(tags))]
    ax_traj = [[plt.subplot(gs_traj[0][i, 0]) for i in range(len(categories))],
               [plt.subplot(gs_traj[1][i, 0]) for i in range(len(categories)-1)]]

    t_ylim        = [0, 1.05]
    t_yticks      = [0, 0.5, 1]
    t_yticklabels = [0, 50, 100]

    t_xlim   = [0, 1850]
    t_xticks = [0, 250, 500, 750, 1000, 1250, 1500, 1750]

    for k in range(len(tags)):

        # load data
        df_traj = pd.read_csv(traj_files[k])
        ids     = np.unique(df_traj['id'])
    
        ## a/b/c.1 -- frequency trajectory
        
        times  = [[] for i in range(len(categories))]
        freqs  = [[] for i in range(len(categories))]
        colors = [[] for i in range(len(categories))]
        
        t_smooth = 100
        for id in ids:
            # temp_times = np.array(df_traj[df_traj['id']==id]['time'])
            # temp_freqs = np.array(df_traj[df_traj['id']==id]['frequency'])

            # x = np.linspace(np.min(temp_times), np.max(temp_times), 500)
            # y = [np.sum(temp_freqs * np.exp(-np.fabs(t - temp_times)/t_smooth) / np.sum(np.exp(-np.fabs(t - temp_times)/t_smooth))) for t in x]

            # times.append(x)
            # freqs.append(y)

            idx = tag2catidx[id.split('.')[0]]

            times[idx].append(df_traj[df_traj['id']==id]['time'])
            freqs[idx].append(df_traj[df_traj['id']==id]['frequency'])
            colors[idx].append(catidx2color[idx])
        
        lineprops = {'lw': SIZELINE, 'ls': '-', 'alpha': 0.5 }
        
        for i in range(len(categories)):
            if len(freqs[i])==0:
                continue

            pprops = { 'xlim':        [lim for lim in t_xlim],
                       'xticks':      [xt for xt in t_xticks],
                       'ylim':        [lim for lim in t_ylim],
                       'yticks':      [yt for yt in t_yticks],
                       'yticklabels': [yl for yl in t_yticklabels],
                       'yminorticks': [0.25, 0.75],
                       'xlabel':      'Time (days after infection)',
                       'ylabel':      '%s\nmutation frequency (%%)' % categories[i],
                       'axoffset':    0.1,
                       'theme':       'open' }
            
            if i<len(ax_traj[k])-1:
                # pprops['xticks'] = []
                pprops['xticklabels'] = ['' for xt in t_xticks]
                pprops['xlabel'] = ''
                # pprops['hide'] = ['bottom']

            mp.plot(type='line', ax=ax_traj[k][i], x=times[i], y=freqs[i], colors=colors[i], plotprops=lineprops, **pprops)
        
    ## legends and labels
    
    sublabels = ['a', 'b']
    ppts = ['CH505', 'CH848']
    dx = -0.08
    dy =  0.0
    for k in range(len(sublabels)):
        ax_traj[k][0].text(box_traj[k]['left'] + dx, box_traj[k]['top'] + dy, sublabels[k].upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
        # ax_traj[k][0].text(box_traj[k]['left'], box_traj[k]['top'], ppts[k], transform=fig.transFigure, **DEF_LABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, filename, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)


def plot_trajectory_selection_shiv(s_files, traj_files, tags, rms, filename='fig-s-vs-time-shiv'):
    """ Plot inferred selection coefficients and times they were first observed, together with selected trajectories. """
    
    # PLOT FIGURE
    
    ## set up figure grid
    
    w     = SINGLE_COLUMN #SLIDE_WIDTH
    goldh = w / 1.35
    fig   = plt.figure(figsize=(w, goldh))

    left_bd   = 0.12
    right_bd  = 0.95
    x_space   = 0.15
    box_x     = (right_bd - left_bd - 1*x_space) * (4/5)
    top_bd    = 0.95
    bottom_bd = 0.15
    y_space   = 0.20
    box_y     = (top_bd - bottom_bd - 1*y_space)/2
    
    box_s    = [dict(left=left_bd + box_x + x_space, right=right_bd, top=top_bd - k*box_y - k*y_space, bottom=top_bd - (k+1)*box_y - k*y_space) for k in range(2)]
    box_traj = [dict(left=left_bd, right=left_bd + box_x, top=top_bd - k*box_y - k*y_space, bottom=top_bd - (k+1)*box_y - k*y_space) for k in range(2)]
    
    ## a -- selection coefficients and times that they were first observed
    
    gs_s = [gridspec.GridSpec(1, 1, **box_s[k]) for k in range(len(tags))]
    ax_s = [plt.subplot(gs_s[k][0, 0]) for k in range(len(tags))]

    gs_traj = [gridspec.GridSpec(1, 1, **box_traj[k]) for k in range(len(tags))]
    ax_traj = [plt.subplot(gs_traj[k][0, 0]) for k in range(len(tags))]
    
    s_ylim        = [-0.01, 0.03]
    s_yticks      = [-0.01, 0, 0.01, 0.02, 0.03]
    s_yticklabels = [int(i) for i in 100*np.array(s_yticks)]

    s_xlim = [0, 4]

    t_ylim        = [0, 1.05]
    t_yticks      = [0, 0.5, 1]
    t_yticklabels = [0, 50, 100]

    t_xlim   = [[0, 450], [0, 900]] #[0, 700]]
    t_xticks = [[0, 150, 300, 450], [0, 300, 600, 900]] #[0, 150, 300, 450, 600, 700]]

    ab_palette  = sns.color_palette(palette='flare', n_colors=5)
    alt_palette = sns.husl_palette(3)
    
    for k in range(len(tags)):

        # load data
        df_s = pd.read_csv(s_files[k])
        df_traj = pd.read_csv(traj_files[k])

        # joint RMs
        # muts = np.unique(df_s['mutation'])
        # for m in muts:
        #     df_temp = df_s[df_s['mutation']==m]
        #     if rms[k] in df_temp['individual'].values:
        #         if 'RMs' not in df_temp['individual'].values:
        #             print(df_temp.head())
        #         df_s.loc[df_s['mutation']==m, 'selection'] = float(df_s[(df_s['mutation']==m) & (df_s['individual']=='RMs')].iloc[0]['selection'])

        df_s = df_s[df_s['individual']==rms[k]]
        df_traj = df_traj[df_traj['individual']==rms[k]]

        ids     = np.unique(df_traj['id'])
        cat     = []
        c_color = []
        c_label = []
        if tags[k]=='SHIV.CH505':
            cat     = ['increase_VL', 'bnAbs_resistance', 'glycan']
            c_color = [C_MPL, ab_palette[0], alt_palette[2]]
            c_label = ['VL', 'bnAb', 'glycan']
        elif tags[k]=='SHIV.CH848':
            cat     = ['bnAbs_resistance', 'glycan']
            c_color = [ab_palette[0], alt_palette[2]]
            c_label = ['bnAb', 'glycan']

        x_std  = 0.1
        x_vals = []
        s_vals = []
        for c in cat:
            if c=='glycan':
                df_sub = df_s[df_s['types'].str.contains('glycan')]
            else:
                df_sub = df_s[df_s['types']==c]
            s_vals.append(df_sub['selection_RMs'])
            x_vals.append(np.random.normal(cat.index(c)+0.5, x_std, len(s_vals[-1])))
    
        ## a/b/c.1 -- frequency trajectory
        
        times  = []
        freqs  = []
        colors = []
        
        t_smooth = 100
        for id in ids:
            # temp_times = np.array(df_traj[df_traj['id']==id]['time'])
            # temp_freqs = np.array(df_traj[df_traj['id']==id]['frequency'])

            # x = np.linspace(np.min(temp_times), np.max(temp_times), 500)
            # y = [np.sum(temp_freqs * np.exp(-np.fabs(t - temp_times)/t_smooth) / np.sum(np.exp(-np.fabs(t - temp_times)/t_smooth))) for t in x]

            # times.append(x)
            # freqs.append(y)

            df_temp = df_traj[df_traj['id']==id]
            idx = -1
            if df_temp['types'].iloc[0] in cat:
                idx = cat.index(df_temp['types'].iloc[0])

            times.append(df_temp['date'])
            freqs.append(df_temp['frequency'])
            colors.append(c_color[idx])
        
        lineprops = {'lw': SIZELINE, 'ls': '-', 'alpha': 0.5 }
        
        pprops = { 'xlim':        [lim for lim in t_xlim[k]],
                   'xticks':      [xt for xt in t_xticks[k]],
                   'ylim':        [lim for lim in t_ylim],
                   'yticks':      [yt for yt in t_yticks],
                   'yticklabels': [yl for yl in t_yticklabels],
                   'yminorticks': [0.25, 0.75],
                   'xlabel':      'Time (days after infection)',
                   'ylabel':      'Frequency (%)',
                   'axoffset':    0.1,
                   'theme':       'open' }

        mp.plot(type='line', ax=ax_traj[k], x=times, y=freqs, colors=colors, plotprops=lineprops, **pprops)
        
        ## a/b/c.2 -- inferred selection coefficient

        scatterprops = dict(lw=0, s=SMALLSIZEDOT, marker='o', alpha=0.35)
        errorprops   = dict(lw=AXWIDTH, markersize=SMALLSIZEDOT, fmt='.', elinewidth=AXWIDTH, capthick=0, capsize=0, mew=AXWIDTH, alpha=1)

        pprops = { 'xlim':        [lim for lim in s_xlim],
                   'xticks':      np.array(range(len(cat)))+0.5,
                   'xticklabels': c_label,
                   'ylim':        [lim for lim in s_ylim],
                   'yticks':      [yt for yt in s_yticks],
                   'yticklabels': [yl for yl in s_yticklabels],
                #    'yminorticks': yminorticks,
                   'ylabel':      'Inferred selection\ncoefficients (%)',
                   'axoffset':    0.1,
                   'theme':       'open' }

        mp.scatter(ax=ax_s[k], x=x_vals, y=s_vals, colors=c_color, plotprops=scatterprops, **pprops)

        
        mp.plot(type='error', ax=ax_s[k], x=[[i+0.5] for i in range(len(c_color))], y=[[np.mean(sv)] for sv in s_vals], yerr=[[np.std(sv)] for sv in s_vals],
                edgecolor=[BKCOLOR for i in c_color], facecolor=[i for i in c_color], plotprops=errorprops, **pprops)
        
        for tick in ax_s[k].get_xticklabels():
            tick.set_rotation(90)

    ## legends and labels
    
    sublabels = ['a', 'b']
    ppts = ['SHIV.CH505 (%s)' % rms[0], 'SHIV.CH848 (%s)' % rms[1]]
    dx = -0.09
    dy =  0.02
    for k in range(len(sublabels)):
        ax_traj[k].text(box_traj[k]['left'] + dx, box_traj[k]['top'] + dy, sublabels[k].upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
        ax_traj[k].text(box_traj[k]['left'], box_traj[k]['top'], ppts[k], transform=fig.transFigure, **DEF_LABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, filename, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)


def plot_selection_vs_rms(s_files, tags, filename='fig-s-vs-rms'):
    """ Plot inferred selection coefficients and times they were first observed, together with selected trajectories. """
    
    # PLOT FIGURE
    
    ## set up figure grid
    
    w     = DOUBLE_COLUMN #SLIDE_WIDTH
    goldh = w / 2.5
    fig   = plt.figure(figsize=(w, goldh))

    left_bd   = 0.08
    right_bd  = 0.95
    x_space   = 0.15
    box_x     = (right_bd - left_bd - 1*x_space)/2
    top_bd    = 0.95
    bottom_bd = 0.15
    
    box_s = [dict(left=left_bd + k*box_x + k*x_space, right=left_bd + (k+1)*box_x + k*x_space, top=top_bd, bottom=bottom_bd) for k in range(2)]
    
    ## a -- selection coefficients
    
    gs_s = [gridspec.GridSpec(1, 1, **box_s[k]) for k in range(len(tags))]
    ax_s = [plt.subplot(gs_s[k][0, 0]) for k in range(len(tags))]
    
    s_ylim        = [-1, 3]
    s_yticks      = [-1, 0, 1, 2, 3]
    s_yticklabels = [int(i) for i in np.array(s_yticks)]

    s_xlim = [0, 7]
    
    for k in range(len(tags)):

        # load data
        df_s = pd.read_csv(s_files[k])

        n_rms  = np.sort(np.unique(df_s['num_RMs']))
        s_vals = []
        x_vals = []
        x_std  = 0.1
        for n in n_rms:
            s_vals.append(df_s[df_s['num_RMs']==n]['selection'])
            x_vals.append(np.random.normal(n-0.5, x_std, len(s_vals[-1])))
        
        ## a/b/c.2 -- inferred selection coefficients

        scatterprops = dict(lw=0, s=SMALLSIZEDOT, marker='o', alpha=0.25)
        errorprops   = dict(lw=AXWIDTH, markersize=SMALLSIZEDOT, fmt='.', elinewidth=AXWIDTH, capthick=0, capsize=0, mew=AXWIDTH, alpha=1)

        pprops = { 'xlim':        [lim for lim in s_xlim],
                   'xticks':      np.array(range(len(n_rms)))+0.5,
                   'xticklabels': n_rms,
                   'ylim':        [lim for lim in s_ylim],
                   'yticks':      [yt for yt in s_yticks],
                   'yticklabels': [yl for yl in s_yticklabels],
                #    'yminorticks': yminorticks,
                   'xlabel':      'Number of RMs in which mutation is observed',
                   'ylabel':      'Inferred selection coefficients (%)',
                   'axoffset':    0.1,
                   'theme':       'open' }

        mp.scatter(ax=ax_s[k], x=x_vals, y=s_vals, colors=[LCOLOR for i in range(len(n_rms))], plotprops=scatterprops, **pprops)

        mp.plot(type='error', ax=ax_s[k], x=[[i+0.5] for i in range(len(n_rms))], y=[[np.mean(sv)] for sv in s_vals], yerr=[[np.std(sv)] for sv in s_vals],
                edgecolor=[BKCOLOR for i in range(len(n_rms))], facecolor=[BKCOLOR for i in range(len(n_rms))], plotprops=errorprops, **pprops)
        
        dashlineprops = {'lw': SIZELINE * 2.0, 'ls': ':', 'alpha': 0.5, 'color': BKCOLOR }
        ax_s[k].axhline(y=0, **dashlineprops)
        
        # for tick in ax_s[k].get_xticklabels():
        #     tick.set_rotation(90)

    ## legends and labels
    
    sublabels = ['a', 'b']
    ppts = ['CH505', 'CH848']
    dx = -0.04
    dy =  0.02
    for k in range(len(sublabels)):
        ax_s[k].text(box_s[k]['left'] + dx, box_s[k]['top'] + dy, sublabels[k].upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
        # ax_s[k].text(box_s[k]['left'], box_s[k]['top'], ppts[k], transform=fig.transFigure, **DEF_LABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, filename, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)


def plot_fitness_comparison(f_files, tags, filename='fig-f-compare'):
    """ Plot inferred selection coefficients and times they were first observed, together with selected trajectories. """
    
    # PLOT FIGURE
    
    ## set up figure grid
    
    w     = SINGLE_COLUMN #SLIDE_WIDTH
    goldh = w * 1.6
    fig   = plt.figure(figsize=(w, goldh))

    left_bd   = 0.15
    top_bd    = 0.98
    bottom_bd = 0.08
    y_space   = 0.08
    box_y     = (top_bd - bottom_bd - 1*y_space)/2
    right_bd  = left_bd + (box_y * (goldh / w))

    print(right_bd)
    
    box_f = [dict(left=left_bd, right=right_bd, top=top_bd - k*box_y - k*y_space, bottom=top_bd - (k+1)*box_y - k*y_space) for k in range(2)]
    
    ## a -- selection coefficients and times that they were first observed
    
    gs_f = [gridspec.GridSpec(1, 1, **box_f[k]) for k in range(len(tags))]
    ax_f = [plt.subplot(gs_f[k][0, 0]) for k in range(len(tags))]

    f_lim    = [-2, 40]
    f_ticks  = [0, 10, 20, 30, 40]
    f_mticks = [5, 15, 25, 35]
    
    for k in range(len(tags)):

        # load data
        df_f = pd.read_csv(f_files[k])

        f_shiv_vals = []
        f_hiv_vals  = []

        df_f_sub = df_f[df_f['breadth']==True]
        f_shiv_vals.append(np.array(df_f_sub['F_jointedRMs'], float))
        f_hiv_vals.append(np.array(df_f_sub['F_%s' % tags[k]], float))

        df_f_sub = df_f[df_f['breadth']==False]
        f_shiv_vals.append(np.array(df_f_sub['F_jointedRMs'], float))
        f_hiv_vals.append(np.array(df_f_sub['F_%s' % tags[k]], float))
    
        ## fitness comparison
        
        scatterprops = dict(lw=0, s=SMALLSIZEDOT, marker='o', alpha=0.4, clip_on=False)
        
        pprops = { 'xlim':        [lim for lim in f_lim],
                   'xticks':      [xt for xt in f_ticks],
                   'xminorticks': [xt for xt in f_mticks],
                   'ylim':        [lim for lim in f_lim],
                   'yticks':      [yt for yt in f_ticks],
                   'yminorticks': [yt for yt in f_mticks],
                   'xlabel':      'SHIV.%s fitness gain (%%)' % tags[k],
                   'ylabel':      '%s fitness gain (%%)' % tags[k],
                   'theme':       'open' }

        mp.plot(type='scatter', ax=ax_f[k], x=f_shiv_vals, y=f_hiv_vals, colors=[C_BNAB, C_NEU], plotprops=scatterprops, **pprops)

        print(tags[k], st.pearsonr(df_f['F_jointedRMs'], df_f['F_%s' % tags[k]]))

        ## legend

        if k==0:
            legendprops = dict(lw=0, s=SMALLSIZEDOT, marker='o', alpha=1, clip_on=False)
            x_legend = 27
            y_legend = 1
            dx_legend = 1
            dy_legend = 2
            label_legend = ['Developed bnAbs', 'No bnAbs']
            color_legend = [C_BNAB, C_NEU]
            for i in range(2):
                mp.scatter(ax=ax_f[k], x=[x_legend], y=[y_legend - i*dy_legend + 0.5], colors=[color_legend[i]], plotprops=legendprops, **pprops)
                ax_f[k].text(x_legend + dx_legend, y_legend - i*dy_legend, label_legend[i], **DEF_LABELPROPS)


    ## legends and labels
    
    sublabels = ['a', 'b']
    ppts = ['CH505', 'CH848']
    dx = -0.08
    dy =  0.0
    for k in range(len(sublabels)):
        ax_f[k].text(box_f[k]['left'] + dx, box_f[k]['top'] + dy, sublabels[k].upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
        # ax_f[k].text(box_f[k]['left'] - 0.02, box_f[k]['top'],      ppts[k],              transform=fig.transFigure, **DEF_LABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, filename, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)


def plot_fitness_comparison_horizontal(f_files, tags, filename='fig-f-compare-h', use_breadth=True):
    """ Plot inferred selection coefficients and times they were first observed, together with selected trajectories. """
    
    # PLOT FIGURE
    
    ## set up figure grid
    
    w     = DOUBLE_COLUMN #SLIDE_WIDTH
    goldh = w / 2.8
    fig   = plt.figure(figsize=(w, goldh))

    top_bd    = 0.95
    bottom_bd = 0.15
    x_space   = 0.10
    box_y     = top_bd - bottom_bd
    x_size    = box_y * (goldh / w)
    x_pad     = 1 - (2*x_size + x_space)
    left_bd   = x_pad/2
    right_bd  = 1 - x_pad/2

    print(x_pad)
    
    box_f = [dict(left=left_bd+k*(x_size+x_space), right=left_bd+k*(x_size+x_space)+x_size, top=top_bd, bottom=bottom_bd) for k in range(2)]
    
    ## a -- selection coefficients and times that they were first observed
    
    gs_f = [gridspec.GridSpec(1, 1, **box_f[k]) for k in range(len(tags))]
    ax_f = [plt.subplot(gs_f[k][0, 0]) for k in range(len(tags))]

    f_lim    = [-2, 40]
    f_ticks  = [0, 10, 20, 30, 40]
    f_mticks = [5, 15, 25, 35]
    
    for k in range(len(tags)):

        # load data
        df_f = pd.read_csv(f_files[k])

        f_shiv_vals = []
        f_hiv_vals  = []
        c_vals      = []

        if use_breadth:
            df_f_sub = df_f[df_f['breadth']==True]
            f_shiv_vals.append(np.array(df_f_sub['F_jointedRMs'], float))
            f_hiv_vals.append(np.array(df_f_sub['F_%s' % tags[k]], float))

            df_f_sub = df_f[df_f['breadth']==False]
            f_shiv_vals.append(np.array(df_f_sub['F_jointedRMs'], float))
            f_hiv_vals.append(np.array(df_f_sub['F_%s' % tags[k]], float))

            c_vals = [C_BNAB, C_NEU]

        else:
            f_shiv_vals.append(np.array(df_f['F_jointedRMs'], float))
            f_hiv_vals.append(np.array(df_f['F_%s' % tags[k]], float))
            c_vals = [C_NEU]
    
        ## fitness comparison
        
        scatterprops = dict(lw=0, s=SMALLSIZEDOT, marker='o', alpha=0.4, clip_on=False)
        
        pprops = { 'xlim':        [lim for lim in f_lim],
                   'xticks':      [xt for xt in f_ticks],
                   'xminorticks': [xt for xt in f_mticks],
                   'ylim':        [lim for lim in f_lim],
                   'yticks':      [yt for yt in f_ticks],
                   'yminorticks': [yt for yt in f_mticks],
                   'xlabel':      'SHIV.%s fitness gain (RM model, %%)' % (tags[k]),
                   'ylabel':      'SHIV.%s fitness gain (human model, %%)' % (tags[k]),
                   'theme':       'open' }

        mp.plot(type='scatter', ax=ax_f[k], x=f_shiv_vals, y=f_hiv_vals, colors=c_vals, plotprops=scatterprops, **pprops)

        print(tags[k], st.pearsonr(df_f['F_jointedRMs'], df_f['F_%s' % tags[k]]))

        ## legend

        if k==0 and use_breadth:
            legendprops = dict(lw=0, s=SMALLSIZEDOT, marker='o', alpha=1, clip_on=False)
            x_legend = 27
            y_legend = 1
            dx_legend = 1
            dy_legend = 2
            label_legend = ['Developed bnAbs', 'No bnAbs']
            color_legend = [C_BNAB, C_NEU]
            for i in range(2):
                mp.scatter(ax=ax_f[k], x=[x_legend], y=[y_legend - i*dy_legend + 0.5], colors=[color_legend[i]], plotprops=legendprops, **pprops)
                ax_f[k].text(x_legend + dx_legend, y_legend - i*dy_legend, label_legend[i], **DEF_LABELPROPS)


    ## legends and labels
    
    sublabels = ['a', 'b']
    ppts = ['CH505', 'CH848']
    dx = -0.05
    dy =  0.0
    for k in range(len(sublabels)):
        ax_f[k].text(box_f[k]['left'] + dx, box_f[k]['top'] + dy, sublabels[k].upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
        # ax_f[k].text(box_f[k]['left'] - 0.02, box_f[k]['top'],      ppts[k],              transform=fig.transFigure, **DEF_LABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, filename, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)


def plot_fitness_gain_v_time(f_files, tags, t_breadth=None, use_breadth=True, filename='fig-gain'):
    """ Plot inferred selection coefficients and times they were first observed, together with selected trajectories. """
    
    # PLOT FIGURE
    
    ## set up figure grid
    
    w     = DOUBLE_COLUMN #SLIDE_WIDTH
    goldh = w / 3.2
    fig   = plt.figure(figsize=(w, goldh))

    top_bd    = 0.95
    bottom_bd = 0.15
    x_space   = 0.10
    left_bd   = 0.10
    right_bd  = 0.95
    x_size    = (right_bd - left_bd - x_space)/2
    
    box_f = [dict(left=left_bd+k*(x_size+x_space), right=left_bd+k*(x_size+x_space)+x_size, top=top_bd, bottom=bottom_bd) for k in range(2)]
    
    ## a -- selection coefficients and times that they were first observed
    
    gs_f = [gridspec.GridSpec(1, 1, **box_f[k]) for k in range(len(tags))]
    ax_f = [plt.subplot(gs_f[k][0, 0]) for k in range(len(tags))]

    f_lim    = [0, 40]
    f_ticks  = [0, 10, 20, 30, 40]
    f_mticks = [5, 15, 25, 35]
    
    x_lim    = [[0, 800], [0, 900]]
    x_ticks  = [[0, 200, 400, 600, 800], [0, 200, 400, 600, 800]]
    x_mticks = [[100, 300, 500, 700], [100, 300, 500, 700, 900]]

    for k in range(len(tags)):

        # load data
        df_f = pd.read_csv(f_files[k])

        t_vals = []
        f_vals = []
        f_stds = []
        c_vals = []
        rm_list = []

        t_broad = []
        f_broad = []
        t_narrow = []
        f_narrow = []

        if use_breadth:
            df_f_sub = df_f[df_f['breadth']==True]
            rms = [rm for rm in np.unique(df_f_sub['RM_id']) if 'RM' in rm]
            for rm in rms:
                t_vals.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['date'], float))
                f_vals.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['mean_F'], float))
                f_stds.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['std_F'], float))
                if rm=='RM6072':
                    c_vals.append(C_MPL)
                else:
                    c_vals.append(C_BNAB)
                    t_broad.append(t_vals[-1])
                    f_broad.append(f_vals[-1])
                rm_list.append(rm)

            df_f_sub = df_f[df_f['breadth']==False]
            rms = [rm for rm in np.unique(df_f_sub['RM_id']) if 'RM' in rm]
            for rm in rms:
                t_vals.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['date'], float))
                f_vals.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['mean_F'], float))
                f_stds.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['std_F'], float))
                c_vals.append(C_NEU)
                rm_list.append(rm)
                t_narrow.append(t_vals[-1])
                f_narrow.append(f_vals[-1])

        else:
            rms = [rm for rm in np.unique(df_f_sub['RM_id']) if 'RM' in rm]
            for rm in rms:
                t_vals.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['date'], float))
                f_vals.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['mean_F'], float))
                f_stds.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['std_F'], float))
                c_vals.append(C_NEU)
                rm_list.append(rm)

        # ## diminishing returns epistasis (power law) curve fitting

        # t_fit = np.arange(0, np.max([np.max(_t) for _t in t_vals]), 1)
        # f_fit_broad = []
        # f_fit_narrow = []

        # if use_breadth:

        #     def df_power_law(t, ab):
        #         return ((ab[1]*t + 1)**ab[0]) - 1
            
        #     def loss(ab, *args):
        #         t  = args[0]
        #         df = args[1]
        #         return df - df_power_law(t, ab)
            
        #     ab_fit = []

        #     print(tags[k])
        #     print('\talpha\t\tbeta')

        #     for label, t_data, f_data in [['broad', t_broad, f_broad], ['narrow', t_narrow, f_narrow]]:
        #         t_data = np.concatenate(tuple(t_data), axis=None)
        #         f_data = np.concatenate(tuple(f_data), axis=None)

        #         ab_data = so.leastsq(loss, (1, 0), args=(t_data, f_data))[0]
        #         ab_fit.append(ab_data)

        #         print('%s\t%.2e\t%.2e' % (label, ab_data[0], ab_data[1]))
            
        #     print('')

        #     f_fit_broad = df_power_law(t_fit, ab_fit[0])
        #     f_fit_narrow = df_power_law(t_fit, ab_fit[1])
        
    
        ## fitness gain
        
        errorprops    = dict(lw=0, markersize=SMALLSIZEDOT, fmt='.', elinewidth=AXWIDTH, capthick=0, capsize=0, mew=0, alpha=1, clip_on=False)
        lineprops     = {'lw': SIZELINE, 'ls': '-', 'alpha': 0.6}
        dashlineprops = {'lw': SIZELINE, 'ls': '--', 'alpha': 0.6}
        
        pprops = { 'xlim':        [lim for lim in x_lim[k]],
                   'xticks':      [xt for xt in x_ticks[k]],
                   'xminorticks': [xt for xt in x_mticks[k]],
                   'ylim':        [lim for lim in f_lim],
                   'yticks':      [yt for yt in f_ticks],
                   'yminorticks': [yt for yt in f_mticks],
                   'xlabel':      'Time (days after infection)',
                   'ylabel':      '%s fitness gain (%%)' % tags[k],
                   'theme':       'open' }

        mp.line(              ax=ax_f[k], x=t_vals, y=f_vals,              colors=c_vals, plotprops=lineprops,  **pprops)
        mp.plot(type='error', ax=ax_f[k], x=t_vals, y=f_vals, yerr=f_stds, colors=c_vals, plotprops=errorprops, **pprops)
        
        # mp.line(ax=ax_f[k], x=[t_fit, t_fit], y=[f_fit_broad, f_fit_narrow], colors=[C_BNAB, C_NEU], plotprops=dashlineprops,  **pprops)

        ## highlight times when breadth established

        scatterprops = dict(lw=AXWIDTH, s=SMALLSIZEDOT*10.0, marker='o', alpha=1)

        if t_breadth is not None:
            for rm in t_breadth:
                if rm in rm_list:
                    rm_idx = rm_list.index(rm)
                    t_idx  = np.argmin(np.fabs(t_breadth[rm] - t_vals[rm_idx]))
                    t_val  = [t_vals[rm_idx][t_idx]]
                    f_val  = [f_vals[rm_idx][t_idx]]
                    mp.scatter(ax=ax_f[k], x=[t_val], y=[f_val], facecolor=['none'], edgecolor=[C_BNAB], plotprops=scatterprops, **pprops)

        ## legend

        if k==0 and use_breadth:
            legendprops = dict(lw=0, s=SMALLSIZEDOT, marker='o', alpha=1, clip_on=False)
            x_legend = 560
            y_legend = 8.5 - 2.5
            dx_legend = 25
            dy_legend = 2.5
            label_legend = ['Developed bnAbs', 'RM6072', 'No bnAbs']
            color_legend = [C_BNAB, C_MPL, C_NEU]
            for i in range(3):
                mp.scatter(ax=ax_f[k], x=[x_legend], y=[y_legend - i*dy_legend + 0.5], colors=[color_legend[i]], plotprops=legendprops, **pprops)
                ax_f[k].text(x_legend + dx_legend, y_legend - i*dy_legend, label_legend[i], **DEF_LABELPROPS)

            if t_breadth is not None:
                mp.scatter(ax=ax_f[k], x=[x_legend], y=[y_legend + dy_legend + 0.5], facecolor=['none'], edgecolor=[C_BNAB], plotprops=scatterprops, **pprops)
                ax_f[k].text(x_legend + dx_legend, y_legend + dy_legend, 'Breadth established', **DEF_LABELPROPS)

            # if use_breadth:
            #     mp.line(ax=ax_f[k], x=[[x_legend-13, x_legend+11]], y=[[y_legend - 3*dy_legend + 0.5, y_legend - 3*dy_legend + 0.5]], colors=[BKCOLOR], plotprops=dashlineprops, **pprops)
            #     ax_f[k].text(x_legend + dx_legend, y_legend - 3*dy_legend, 'Power law fit', **DEF_LABELPROPS)


    ## legends and labels
    
    sublabels = ['a', 'b']
    ppts = ['CH505', 'CH848']
    dx = -0.05
    dy =  0.0
    for k in range(len(sublabels)):
        ax_f[k].text(box_f[k]['left'] + dx, box_f[k]['top'] + dy, sublabels[k].upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
        # ax_f[k].text(box_f[k]['left'] - 0.02, box_f[k]['top'],      ppts[k],              transform=fig.transFigure, **DEF_LABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, filename, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)


def plot_fitness_gain_v_time_categories(f_files, tags, t_breadth=None, use_breadth=True, filename='fig-gain-categories'):
    """ Plot inferred selection coefficients and times they were first observed, together with selected trajectories. """
    
    # PLOT FIGURE
    
    ## set up figure grid
    
    w     = DOUBLE_COLUMN #SLIDE_WIDTH
    goldh = w / 1
    fig   = plt.figure(figsize=(w, goldh))

    top_bd    = 0.97
    bottom_bd = 0.08
    x_space   = 0.10
    left_bd   = 0.10
    right_bd  = 0.95
    x_size    = (right_bd - left_bd - x_space)/2
    
    box_f = [dict(left=left_bd+k*(x_size+x_space), right=left_bd+k*(x_size+x_space)+x_size, top=top_bd, bottom=bottom_bd) for k in range(2)]
    
    ## a -- fitness gain over time by category

    categories = ['Abs-resistance', 'glycan', 'reversion', 'VL-enhancing']
    cat2label  = ['Ab resistance', 'glycan', 'reversion', 'VL enhancing']
    
    gs_f = [gridspec.GridSpec(len(categories), 1, **box_f[k]) for k in range(len(tags))]
    ax_f = [[plt.subplot(gs_f[0][i, 0]) for i in range(len(categories))],
            [plt.subplot(gs_f[1][i, 0]) for i in range(len(categories)-1)]]

    f_lim    = [0, 15]
    f_ticks  = [0, 5, 10, 15]
    f_mticks = [2.5, 7.5, 12.5]
    
    x_lim    = [[0, 800], [0, 900]]
    x_ticks  = [[0, 200, 400, 600, 800], [0, 200, 400, 600, 800]]
    x_mticks = [[100, 300, 500, 700], [100, 300, 500, 700, 900]]

    for k in range(len(tags)):
        # load data
        df_f = pd.read_csv(f_files[k])

        t_vals_c = []
        f_vals_c = []
        f_stds_c = []
        c_vals_c = []
        rm_list_c = []

        t_broad_c = []
        f_broad_c = []
        t_narrow_c = []
        f_narrow_c = []

        for c in categories:
            t_vals = []
            f_vals = []
            f_stds = []
            c_vals = []
            rm_list = []

            t_broad = []
            f_broad = []
            t_narrow = []
            f_narrow = []

            if use_breadth:
                df_f_sub = df_f[(df_f['breadth']==True) & (df_f['Type']==c)]
                rms = [rm for rm in np.unique(df_f_sub['RM_id']) if 'RM' in rm]
                for rm in rms:
                    t_vals.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['date'], float))
                    f_vals.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['mean_F'], float))
                    f_stds.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['std_F'], float))
                    if rm=='RM6072':
                        c_vals.append(C_MPL)
                    else:
                        c_vals.append(C_BNAB)
                        t_broad.append(t_vals[-1])
                        f_broad.append(f_vals[-1])
                    rm_list.append(rm)

                df_f_sub = df_f[(df_f['breadth']==False) & (df_f['Type']==c)]
                rms = [rm for rm in np.unique(df_f_sub['RM_id']) if 'RM' in rm]
                for rm in rms:
                    t_vals.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['date'], float))
                    f_vals.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['mean_F'], float))
                    f_stds.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['std_F'], float))
                    if rm=='RM6072':
                        c_vals.append(C_MPL)
                    else:
                        c_vals.append(C_NEU)
                    rm_list.append(rm)
                    t_narrow.append(t_vals[-1])
                    f_narrow.append(f_vals[-1])

            else:
                df_f_sub = df_f[df_f['Type']==c]
                rms = [rm for rm in np.unique(df_f_sub['RM_id']) if 'RM' in rm]
                for rm in rms:
                    t_vals.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['date'], float))
                    f_vals.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['mean_F'], float))
                    f_stds.append(np.array(df_f_sub[df_f_sub['RM_id']==rm]['std_F'], float))
                    c_vals.append(C_NEU)
                    rm_list.append(rm)

            t_vals_c.append(t_vals)
            f_vals_c.append(f_vals)
            f_stds_c.append(f_stds)
            c_vals_c.append(c_vals)

            t_broad_c.append(t_broad)
            f_broad_c.append(f_broad)

        # ## diminishing returns epistasis (power law) curve fitting

        # t_fit = np.arange(0, np.max([np.max(_t) for _t in t_vals]), 1)
        # f_fit_broad = []
        # f_fit_narrow = []

        # if use_breadth:

        #     def df_power_law(t, ab):
        #         return ((ab[1]*t + 1)**ab[0]) - 1
            
        #     def loss(ab, *args):
        #         t  = args[0]
        #         df = args[1]
        #         return df - df_power_law(t, ab)
            
        #     ab_fit = []

        #     print(tags[k])
        #     print('\talpha\t\tbeta')

        #     for label, t_data, f_data in [['broad', t_broad, f_broad], ['narrow', t_narrow, f_narrow]]:
        #         t_data = np.concatenate(tuple(t_data), axis=None)
        #         f_data = np.concatenate(tuple(f_data), axis=None)

        #         ab_data = so.leastsq(loss, (1, 0), args=(t_data, f_data))[0]
        #         ab_fit.append(ab_data)

        #         print('%s\t%.2e\t%.2e' % (label, ab_data[0], ab_data[1]))
            
        #     print('')

        #     f_fit_broad = df_power_law(t_fit, ab_fit[0])
        #     f_fit_narrow = df_power_law(t_fit, ab_fit[1])
        
    
        ## fitness gain
        
        errorprops    = dict(lw=0, markersize=SMALLSIZEDOT, fmt='.', elinewidth=AXWIDTH, capthick=0, capsize=0, mew=0, alpha=1, clip_on=False)
        lineprops     = {'lw': SIZELINE, 'ls': '-', 'alpha': 0.6}
        dashlineprops = {'lw': SIZELINE, 'ls': '--', 'alpha': 0.6}

        for i in range(len(categories)):

            if len(t_vals_c[i]) == 0:
                continue

            pprops = { 'xlim':        [lim for lim in x_lim[k]],
                       'xticks':      [xt for xt in x_ticks[k]],
                       'xminorticks': [xt for xt in x_mticks[k]],
                       'ylim':        [lim for lim in f_lim],
                       'yticks':      [yt for yt in f_ticks],
                       'yminorticks': [yt for yt in f_mticks],
                       'xlabel':      'Time (days after infection)',
                       'ylabel':      'Fitness gain from\n%s mutations (%%)' % cat2label[i],
                       'theme':       'open' }
            
            if i<len(ax_f[k])-1:
                # pprops['xticks'] = []
                pprops['xticklabels'] = ['' for _ in pprops['xticks']]
                # pprops['xminorticks'] = []
                pprops['xlabel'] = ''
                # pprops['hide'] = ['bottom']

            mp.line(              ax=ax_f[k][i], x=t_vals_c[i], y=f_vals_c[i],                   colors=c_vals_c[i], plotprops=lineprops,  **pprops)
            mp.plot(type='error', ax=ax_f[k][i], x=t_vals_c[i], y=f_vals_c[i], yerr=f_stds_c[i], colors=c_vals_c[i], plotprops=errorprops, **pprops)
            
            # mp.line(ax=ax_f[k], x=[t_fit, t_fit], y=[f_fit_broad, f_fit_narrow], colors=[C_BNAB, C_NEU], plotprops=dashlineprops,  **pprops)

        # ## highlight times when breadth established

        # scatterprops = dict(lw=AXWIDTH, s=SMALLSIZEDOT*10.0, marker='o', alpha=1)

        # if t_breadth is not None:
        #     for rm in t_breadth:
        #         if rm in rm_list:
        #             rm_idx = rm_list.index(rm)
        #             t_idx  = np.argmin(np.fabs(t_breadth[rm] - t_vals[rm_idx]))
        #             t_val  = [t_vals[rm_idx][t_idx]]
        #             f_val  = [f_vals[rm_idx][t_idx]]
        #             mp.scatter(ax=ax_f[k], x=[t_val], y=[f_val], facecolor=['none'], edgecolor=[C_BNAB], plotprops=scatterprops, **pprops)

        ## legend

        if k==0 and use_breadth:
            legendprops = dict(lw=0, s=SMALLSIZEDOT, marker='o', alpha=1, clip_on=False)
            x_legend = 580 # 20
            y_legend = 3.0 # 9.5 
            dx_legend = 25/1.5
            dy_legend = 2.5/3 * (15/10)
            label_legend = ['Developed bnAbs', 'RM6072', 'No bnAbs']
            color_legend = [C_BNAB, C_MPL, C_NEU]
            for i in range(3):
                mp.scatter(ax=ax_f[k][0], x=[x_legend], y=[y_legend - i*dy_legend + 0.25], colors=[color_legend[i]], plotprops=legendprops, **pprops)
                ax_f[k][0].text(x_legend + dx_legend, y_legend - i*dy_legend, label_legend[i], **DEF_LABELPROPS)

            # if t_breadth is not None:
            #     mp.scatter(ax=ax_f[k], x=[x_legend], y=[y_legend + dy_legend + 0.5], facecolor=['none'], edgecolor=[C_BNAB], plotprops=scatterprops, **pprops)
            #     ax_f[k].text(x_legend + dx_legend, y_legend + dy_legend, 'Breadth established', **DEF_LABELPROPS)

            # if use_breadth:
            #     mp.line(ax=ax_f[k][0], x=[[x_legend-13, x_legend+11]], y=[[y_legend - 3*dy_legend + 0.5, y_legend - 3*dy_legend + 0.5]], colors=[BKCOLOR], plotprops=dashlineprops, **pprops)
            #     ax_f[k][0].text(x_legend + dx_legend, y_legend - 3*dy_legend, 'Power law fit', **DEF_LABELPROPS)


    ## legends and labels
    
    sublabels = ['a', 'b']
    ppts = ['CH505', 'CH848']
    dx = -0.05
    dy =  0.0
    for k in range(len(sublabels)):
        ax_f[k][0].text(box_f[k]['left'] + dx, box_f[k]['top'] + dy, sublabels[k].upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
        # ax_f[k].text(box_f[k]['left'] - 0.02, box_f[k]['top'],      ppts[k],              transform=fig.transFigure, **DEF_LABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, filename, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)


#########################
# SUPPLEMENTARY FIGURES #
#########################

def plot_supplementary_figure_example_mpl(**pdata):
    """
    Example evolutionary trajectory for a 50-site system and inferred selection coefficients,
    together with aggregate properties across sampling levels.
    """
    
    # unpack passed data

    n_gen   = pdata['n_gen']
    dg      = pdata['dg']
    N       = pdata['N']
    xfile   = pdata['xfile']
    method  = pdata['method']
    name    = pdata['name']
    hist_ns = pdata['hist_ns']
    hist_dt = pdata['hist_dt']

    n_ben = pdata['n_ben']
    n_neu = pdata['n_neu']
    n_del = pdata['n_del']
    s_ben = pdata['s_ben']
    s_neu = pdata['s_neu']
    s_del = pdata['s_del']

    # load and process data files

    data  = np.loadtxt('%s/data/%s.dat' % (WFS_DIR, xfile))
    times = np.unique(data.T[0])
    x     = []
    for i in range(0, n_gen, dg):
        idx    = data.T[0]==i
        t_data = data[idx].T[2:].T
        t_num  = data[idx].T[1].T
        t_freq = np.einsum('i,ij->j', t_num, t_data) / float(np.sum(t_num))
        x.append(t_freq)
    x = np.array(x).T

    s_true = [s_ben for i in range(n_ben)] + [0 for i in range(n_neu)] + [s_del for i in range(n_del)]
    s_inf  = np.loadtxt('%s/%s_%s.dat' % (SIM_MPL_DIR, xfile.split('wfsim_')[1], method))
    cov    = np.loadtxt('%s/covariance-%s.dat' % (SIM_MPL_DIR, xfile.split('wfsim_')[1]))
    ds     = np.linalg.inv(cov) / N

    # PLOT FIGURE

    ## set up figure grid

    w     = DOUBLE_COLUMN
    goldh = w / 1.08
    fig   = plt.figure(figsize=(w, goldh))

    n_rows = [   n_ben,    n_neu,       n_del]
    offset = [       0,    n_ben, n_ben+n_neu]
    tag    = [   'ben',    'neu',       'del']
    colors = [   C_BEN,    C_NEU,       C_DEL]
    fc     = [C_BEN_LT, C_NEU_LT,    C_DEL_LT]

    htot = 0.80 + 0.01
    dh   = (htot - 0.08) / float(n_ben + n_neu + n_del + 2)

    boxl = [0.06, 0.06, 0.06, 0.06]
    boxr = [0.34, 0.34, 0.34, 0.34]
    boxb = [htot - (dh * n_ben), htot - (dh * (n_ben + n_neu + 1)), htot - (dh * (n_ben + n_neu + n_del + 2))]
    boxt = [               htot,         htot - (dh * (n_ben + 1)), htot - (dh * (n_ben + n_neu + 2))]
    gs   = [0 for k in range(len(tag))]

    ## a -- all trajectories together

    dx = 0.03
    dy = 0.06
    tempgs = gridspec.GridSpec(1, 1)
    tempgs.update(left=boxl[0], right=boxr[0], bottom=htot+dy+0.02, top=0.97)
    ax     = plt.subplot(tempgs[0, 0])

    lineprops = { 'lw' : SIZELINE*1.5, 'ls' : '-', 'alpha' : 0.6 }

    pprops = { 'xticks'      : [0, 50, 100, 150, 200, 250, 300, 350, 400],
               'yticks'      : [0, 1],
               'yminorticks' : [0.25, 0.5, 0.75],
               'nudgey'      : 1.1,
               'xlabel'      : 'Generation',
               'ylabel'      : 'Allele\nfrequency, '+r'$x$',
               'plotprops'   : lineprops,
               'axoffset'    : 0.1,
               'theme'       : 'open' }

    xdat = [range(0, n_gen, dg) for k in range(len(x))]
    ydat = [k for k in x]
    mp.plot(type='line', ax=ax, x=xdat, y=ydat, colors=[LCOLOR for k in range(len(x))], **pprops)

    ## b -- individual beneficial/neutral/deleterious trajectories and selection coefficients

    idx = 0
    for k in range(len(tag)):
        
        ### trajectories
        
        gs[k] = gridspec.GridSpec(n_rows[k], 2)
        gs[k].update(left=boxl[k], right=boxr[k], bottom=boxb[k], top=boxt[k], wspace=0.05)
        ax = [[plt.subplot(gs[k][i, 0]), plt.subplot(gs[k][i, 1])] for i in range(n_rows[k])]
        
        legendprops = {'loc' : 4, 'frameon' : False, 'scatterpoints' : 1, 'handletextpad' : 0.1,
                       'prop' : {'size' : SIZELABEL}, 'ncol' : 1}
        lineprops   = {'lw' : SIZELINE*1.5, 'linestyle' : '-', 'alpha' : 1.0}
        fillprops   = {'lw' : 0, 'alpha' : 0.3, 'interpolate' : True}
        
        pprops = { 'xticks' : [],
                   'yticks' : [],
                   'hide'   : ['top','bottom','left','right'],
                   'theme'  : 'open' }
        
        for i in range(n_rows[k]):
            pprops['colors'] = [colors[k]]
            pprops['xlim']   = [    0,  400]
            pprops['ylim']   = [-0.08, 1.08]
            ydat             = x[offset[k]+i]
            if (i==n_rows[k]-1) and (k==len(tag)-1):
                pprops['xticks']   = [0, 100, 200, 300, 400]
                pprops['xlabel']   = 'Generation'
                pprops['hide']     = ['top','left','right']
                pprops['axoffset'] = 0.3
            mp.line(             ax=ax[i][0], x=[range(0, n_gen, dg)], y=[ydat], plotprops=lineprops, **pprops)
            mp.plot(type='fill', ax=ax[i][0], x=[range(0, n_gen, dg)], y=[ydat], plotprops=fillprops, **pprops)
        
        ### selection coefficient estimates
        
        sprops = {'lw': 0, 's': 9., 'marker': 'o'}
        
        pprops = {'yticks': [],
                  'xticks': [],
                  'hide':   ['top','bottom','left','right'],
                  'theme':  'open' }
        
        for i in range(n_rows[k]):
            pprops['xlim'] = [-0.04, 0.04]
            pprops['ylim'] = [  0.5,  1.5]
            ydat           = [1]
            xdat           = [s_inf[offset[k]+i]]
            xerr           = np.sqrt(ds[offset[k]+i][offset[k]+i])
            if (i==n_rows[k]-1) and (k==len(tag)-1):
                pprops['xticks']      = [-0.04,   -0.02,      0,   0.02,   0.04]
                pprops['xticklabels'] = [  ' ', r'$-2$', r'$0$', r'$2$', r'$4$']
                pprops['xlabel']      = 'Inferred selection\ncoefficient, ' + r'$\hat{s}$' + ' (%)'
                pprops['hide']        = ['top','left','right']
                pprops['axoffset']    = 0.3
            mp.plot(type='error', ax=ax[i][1], x=[xdat], y=[ydat], xerr=[xerr], colors=[colors[k]], **pprops)
            ax[i][1].axvline(x=s_true[idx], ls=':', lw=SIZELINE, color=BKCOLOR)
            idx += 1

    ### bounding boxes

    ax[0][0].text(boxl[0]-0.03,    0.98, 'a'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax[0][0].text(boxl[0]-0.03, boxt[0], 'b'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    lineprops = { 'lw' : AXWIDTH/2., 'ls' : '-', 'alpha' : 1.0 }
    pprops    = { 'xlim' : [0, 1], 'xticks' : [], 'ylim' : [0, 1], 'yticks' : [],
        'hide' : ['top','bottom','left','right'], 'plotprops' : lineprops }
    txtprops = {'ha' : 'right', 'va' : 'center', 'color' : BKCOLOR, 'family' : FONTFAMILY,
        'size' : SIZELABEL, 'rotation' : 90, 'transform' : fig.transFigure}
    ax[0][0].text(boxl[0]-0.01, (boxb[0]+boxt[0])/2.,  'Beneficial', **txtprops)
    ax[0][0].text(boxl[0]-0.01, (boxb[1]+boxt[1])/2.,     'Neutral', **txtprops)
    ax[0][0].text(boxl[0]-0.01, (boxb[2]+boxt[2])/2., 'Deleterious', **txtprops)

    boxprops = {'ec' : BKCOLOR, 'lw' : SIZELINE/2., 'fc' : 'none', 'clip_on' : False, 'zorder' : -100}

    dx  = 0.005
    dxl = 0.000
    dxr = 0.001
    dy  = 0.003
    for k in range(len(tag)):
        ll = boxl[k] + dxl                     # left box left
        lb = boxb[k] - dy                      # left box bottom
        rl = (boxl[k] + boxr[k])/2. + dx + dxl # right box left
        wd = (boxr[k] - boxl[k])/2. - dxr      # box width
        ht = (boxt[k] - boxb[k]) + (2. * dy)   # box height
        
        recL = matplotlib.patches.Rectangle(xy=(ll, lb), width=wd, height=ht, transform=fig.transFigure, **boxprops)
        recL = ax[0][0].add_patch(recL)
        recR = matplotlib.patches.Rectangle(xy=(rl, lb), width=wd, height=ht, transform=fig.transFigure, **boxprops)
        recR = ax[0][0].add_patch(recR)

    ## c -- histogram of selection coefficients with perfect sampling

    ### set up grid

    box_l  = dict(left=0.44, right=0.58, bottom=0.05, top=0.96)
    box_ru = dict(left=0.69, right=0.97, bottom=0.63, top=0.96)
    box_rl = dict(left=0.69, right=0.97, bottom=0.14, top=0.47)

    gs_l  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_l)
    gs_ru = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_ru)
    gs_rl = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_rl)
    ax_l  = plt.subplot(gs_l[ 0, 0])
    ax_ru = plt.subplot(gs_ru[0, 0])
    ax_rl = plt.subplot(gs_rl[0, 0])

    ### plot histogram

    df   = pd.read_csv('%s/MPL_%s_collected_extended.csv.gz' % (SIM_DIR, name), memory_map=True)
    df   = df[df.method==method]
    df_s = df[(df.deltat==hist_dt) & (df.ns==hist_ns)]

    ben_cols = ['s%d' % i for i in range(n_ben)]
    neu_cols = ['s%d' % i for i in range(n_ben, n_ben+n_neu)]
    del_cols = ['s%d' % i for i in range(n_ben+n_neu, n_ben+n_neu+n_del)]

    colors     = [C_BEN, C_NEU, C_DEL]
    tags       = ['beneficial', 'neutral', 'deleterious']
    cols       = [ben_cols, neu_cols, del_cols]
    s_true_loc = [s_ben, s_neu, s_del]

    dashlineprops = { 'lw' : SIZELINE * 2.0, 'ls' : ':', 'alpha' : 0.5, 'color' : BKCOLOR }
    histprops = dict(histtype='bar', lw=SIZELINE/2, rwidth=0.8, ls='solid', alpha=0.7, edgecolor='none',
                     orientation='horizontal')
    pprops = { 'xlim'        : [-0.04, 0.04],
               'xticks'      : [  -0.04,   -0.03,   -0.02,   -0.01,     0.,   0.01,   0.02,   0.03,   0.04],
               'yticklabels' : [r'$-4$', r'$-3$', r'$-2$', r'$-1$', r'$0$', r'$1$', r'$2$', r'$3$', r'$4$'],
               'ylim'        : [0., 0.10],
               'yticks'      : [0., 0.05, 0.10],
               'xlabel'      : 'Frequency',
               'ylabel'      : 'Inferred selection coefficient, ' + r'$\hat{s}$' + ' (%)',
               'bins'        : np.arange(-0.04, 0.04+0.001, 0.001),
               'combine'     : True,
               'plotprops'   : histprops,
               'axoffset'    : 0.1,
               'theme'       : 'boxed' }

    for i in range(len(tags)):
        x = [np.array(df_s[cols[i]]).flatten()]
        tprops = dict(ha='left', va='center', family=FONTFAMILY, size=SIZELABEL, rotation=270, clip_on=False)
        ax_l.text(0.102, s_true_loc[i], r'$s_{%s}$' % (tags[i]), color=colors[i], **tprops)
        ax_l.axhline(y=s_true_loc[i], **dashlineprops)
        if i<len(tags)-1: mp.hist(             ax=ax_l, x=x, colors=[colors[i]], **pprops)
        else:             mp.plot(type='hist', ax=ax_l, x=x, colors=[colors[i]], **pprops)

    ## d, e -- AUCs for inferring beneficial/deleterious mutations

    ns_vals = [10, 20, 30, 40, 50, 80, 100]
    dt_vals = [1, 5, 10, 20, 50]

    AUC_matrix_bn = np.zeros((len(dt_vals), len(ns_vals)))
    AUC_matrix_nd = np.zeros((len(dt_vals), len(ns_vals)))

    for i in range(len(dt_vals)):
        for j in range(len(ns_vals)):
            df_AUC = df[(df.deltat==dt_vals[i]) & (df.ns==ns_vals[j])]
            AUC_matrix_bn[i, j] = np.mean(df_AUC.AUROC_ben)
            AUC_matrix_nd[i, j] = np.mean(df_AUC.AUROC_del)

    pprops = { 'xlim'        : [0, len(dt_vals)],
               'xticks'      : np.arange(len(dt_vals))+0.5,
               'xticklabels' : [int(k) for k in dt_vals],
               'ylim'        : [0, len(ns_vals)],
               'yticks'      : np.arange(len(ns_vals))+0.5,
               'yticklabels' : [int(k) for k in ns_vals],
               'xlabel'      : 'Time between samples, '+r'$\Delta t$' + ' (generations)',
               'ylabel'      : 'Number of sequences per time point, '+r'$n_s$',
               'theme'       : 'boxed' }
    tprops = dict(ha='center', va='center', family=FONTFAMILY, size=SIZELABEL, clip_on=False)

    ax_ru.pcolor(AUC_matrix_bn.T, vmin=0.75, vmax=1.0, cmap='GnBu', alpha=0.75)
    for i in range(len(AUC_matrix_bn)):
        for j in range(len(AUC_matrix_bn[0])):
            tc = 'k'
            if AUC_matrix_bn[i,j]>0.96: tc = 'white'
            ax_ru.text(i+0.5, j+0.5, '%.2f' % (AUC_matrix_bn[i,j]), color=tc, **tprops)
    mp.plot(type='scatter', ax=ax_ru, x=[[-1]], y=[[-1]], colors=[BKCOLOR], **pprops)

    ax_rl.pcolor(AUC_matrix_nd.T, vmin=0.75, vmax=1.0, cmap='GnBu', alpha=0.75)
    for i in range(len(AUC_matrix_nd)):
        for j in range(len(AUC_matrix_nd[0])):
            tc = 'k'
            if AUC_matrix_nd[i,j]>0.96: tc = 'white'
            ax_rl.text(i+0.5, j+0.5, '%.2f' % (AUC_matrix_nd[i,j]), color=tc, **tprops)
    mp.plot(type='scatter', ax=ax_rl, x=[[-1]], y=[[-1]], colors=[BKCOLOR], **pprops)

    ## outside text labels

    tprops = dict(color=BKCOLOR, ha='center', va='center', family=FONTFAMILY, size=SIZELABEL,
                  clip_on=False, transform=fig.transFigure)
    dx = -0.03
    dy =  0.02

    ax_ru.text((box_ru['right']-box_ru['left'])/2+box_ru['left'],  box_l['top']+dy, 'Mean AUROC (beneficial)',  **tprops)
    ax_ru.text((box_rl['right']-box_rl['left'])/2+box_rl['left'], box_rl['top']+dy, 'Mean AUROC (deleterious)', **tprops)

    ax_l.text(  box_l['left']+dx,  box_l['top']+dy, 'c'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_ru.text(box_ru['left']+dx, box_ru['top']+dy, 'd'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_rl.text(box_rl['left']+dx, box_rl['top']+dy, 'e'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/ed-fig-1-%s-mpl%s' % (FIG_DIR, name, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('MPL supplementary example done.')
    

def plot_supplementary_figure_performance(**pdata):
    """
    Comparisons versus existing methods of selection inference.
    """
    
    # unpack data

    test_sets = pdata['test_sets']
    traj_file = pdata['traj_file']
    t_ticks   = pdata['t_ticks']
    n_ben     = pdata['n_ben']
    n_neu     = pdata['n_neu']
    n_del     = pdata['n_del']
    x_ben     = pdata['x_ben']
    y_ben     = pdata['y_ben']
    x_del     = pdata['x_del']
    y_del     = pdata['y_del']
    x_err     = pdata['x_err']
    y_err     = pdata['y_err']
    x_t       = pdata['x_t']
    y_t       = pdata['y_t']

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 1
    goldh   = 0.66 * w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    box_top   = 0.94
    box_left  = 0.15
    box_right = 0.85
    ddy       = 0.11 / hshrink
    dy        = 0.21 / hshrink

    box_cben = dict(left=box_left, right=box_right, bottom=box_top-(1*dy)-(0*ddy), top=box_top-(0*dy)-(0*ddy))
    box_cdel = dict(left=box_left, right=box_right, bottom=box_top-(2*dy)-(1*ddy), top=box_top-(1*dy)-(1*ddy))
    box_rmse = dict(left=box_left, right=box_right, bottom=box_top-(3*dy)-(2*ddy), top=box_top-(2*dy)-(2*ddy))
    
    wspace  = 0.10

    gs_cben = gridspec.GridSpec(1, len(test_sets), wspace=0.10, **box_cben)
    gs_cdel = gridspec.GridSpec(1, len(test_sets), wspace=0.10, **box_cdel)
    gs_rmse = gridspec.GridSpec(1, len(test_sets), wspace=0.10, **box_rmse)
    
    ax_cben = [plt.subplot(gs_cben[0, i]) for i in range(len(test_sets))]
    ax_cdel = [plt.subplot(gs_cdel[0, i]) for i in range(len(test_sets))]
    ax_rmse = [plt.subplot(gs_rmse[0, i]) for i in range(len(test_sets))]
    
    ## set colors and methods list

    hc        = C_MPL
    nc        = C_NEU_LT
    hfc       = '#ffcd5e'
    nfc       = '#f0f0f0'
    methods   = ['MPL', 'FIT', 'LLS', 'CLEAR', 'EandR', 'ApproxWF', 'WFABC',  'IM']
    labels    = ['MPL',   '1',   '2',     '3',     '4',        '5',     '6',   '7']
    colorlist = sns.husl_palette(len(methods)-1)
#    colorbg   = [c + [0.2] for c in colorlist]  # Code Ocean
    colorbg   = [c + tuple([0.2]) for c in colorlist]  # GitHub
    
    lineprops     = { 'lw': SIZELINE*1.2, 'linestyle': '-', 'alpha': 1.0, 'drawstyle': 'steps-mid' }
    fillprops     = { 'lw': 0, 'alpha': 0.2, 'interpolate': True, 'step': 'mid' }
    dashlineprops = { 'lw': SIZELINE*1.2, 'ls': ':', 'alpha': 1.0, 'color': BKCOLOR }
    
    tprops = dict(ha='center', va='top', family=FONTFAMILY, size=SIZELABEL, clip_on=False)
    
    N_BINS = 25

    ## a -- classification of beneficial mutants

    for k in range(len(test_sets)):
        xmin, xmax = -0.14, 0.50
        ymin, ymax =  0, 40
        
        dmax = xmax - xmin
        dval = [np.array(y_ben[k][0]) - np.array(y_ben[k][i]) for i in range(1, len(methods))]
        
        n_bins = N_BINS
        bin_d  = dmax/n_bins
        bins   = np.arange(xmin, xmax+bin_d, bin_d)
        hist_y = []

        for kk in range(len(methods)-1):
            temp_hist_y = [np.sum((bins[i]<=dval[kk]) & (dval[kk]<bins[i+1])) for i in range(len(bins)-1)]
            temp_hist_y.append(np.sum(dval[kk]>=bins[-1]))
            #print(np.sum(temp_hist_y))
            hist_y.append(np.insert(np.append(temp_hist_y, 0), 0, 0))
                    
        bins       = np.insert(np.append(bins, xmax+bin_d), 0, 0)
        xmax      += bin_d
        hist_props = dict(lw=AXWIDTH/2, width=0.9*dmax/n_bins, align='center', orientation='vertical', edgecolor=[BKCOLOR], alpha=0.6)
    
        ax_cben[k].text(xmin/2, 35, 'alternative\nmore\naccurate', color=BKCOLOR, rotation=0, **tprops)
        ax_cben[k].text(xmax/2, 35, 'MPL more accurate', color=BKCOLOR, rotation=0, **tprops)
        ax_cben[k].axvline(x=0, **dashlineprops)

        pprops = { 'colors':      colorlist,
                   'xlim':        [xmin, xmax],
                   'ylim':        [ymin, ymax],
                   'xticks':      [-0.10, 0, 0.10, 0.20, 0.30, 0.40, 0.50],
                   'xticklabels': ['-0.1', '0', '0.1', '0.2', '0.3', '0.4', r'$\geq 0.5$'],
                   'yticks':      [],
                   'xaxstart':    xmin,
                   'xaxend':      xmax,
                   'axoffset':    0,
                   'xlabel':      'Improvement in classification of beneficial alleles (AUROC)',
                   'theme':       'open',
                   'hide':        ['left','right'] }

        if k==0:
            pprops['yticks'] = [0, 10, 20, 30, 40]
            pprops['ylabel'] = 'Number of simulations'
            pprops['hide']   = []
        
        mp.line(             ax=ax_cben[k], x=[bins+(bin_d/2) for kk in range(len(methods)-1)], y=hist_y, plotprops=lineprops, **pprops)
        mp.plot(type='fill', ax=ax_cben[k], x=[bins+(bin_d/2) for kk in range(len(methods)-1)], y=hist_y, plotprops=fillprops, **pprops)
        
        #mp.plot(type='bar', ax=ax_cben[k], x=[bins+(bin_d/2) for kk in range(len(methods)-1)], y=hist_y, plotprops=hist_props, **pprops)
        
        ax_cben[k].text(box_left + (1.04*k+0.5)*(box_right - box_left)/len(test_sets), box_cben['top']+0.03,
                        test_sets[k].split('_')[1].capitalize(),
                        ha='center', va='center', transform=fig.transFigure, clip_on=False, **DEF_LABELPROPS)


    ## b -- classification of deleterious mutants

    for k in range(len(test_sets)):
        xmin, xmax = -0.14, 0.50
        ymin, ymax =  0, 40
        
        dmax = xmax - xmin
        dval = [np.array(y_del[k][0]) - np.array(y_del[k][i]) for i in range(1, len(methods))]
        
        n_bins = N_BINS
        bin_d  = dmax/n_bins
        bins   = np.arange(xmin, xmax+bin_d, bin_d)
        hist_y = []

        for kk in range(len(methods)-1):
            temp_hist_y = [np.sum((bins[i]<=dval[kk]) & (dval[kk]<bins[i+1])) for i in range(len(bins)-1)]
            temp_hist_y.append(np.sum(dval[kk]>=bins[-1]))
            #print(np.sum(temp_hist_y))
            hist_y.append(np.insert(np.append(temp_hist_y, 0), 0, 0))
                    
        bins       = np.insert(np.append(bins, xmax+bin_d), 0, 0)
        xmax      += bin_d
        hist_props = dict(lw=AXWIDTH/2, width=0.9*dmax/n_bins, align='center', orientation='vertical', edgecolor=[BKCOLOR], alpha=0.6)
            
        ax_cdel[k].text(xmin/2, 35, 'alternative\nmore\naccurate', color=BKCOLOR, rotation=0, **tprops)
        ax_cdel[k].text(xmax/2, 35, 'MPL more accurate', color=BKCOLOR, rotation=0, **tprops)
        ax_cdel[k].axvline(x=0, **dashlineprops)

        pprops = { 'colors':      colorlist,
                   'xlim':        [xmin, xmax],
                   'ylim':        [ymin, ymax],
                   'xticks':      [-0.10, 0, 0.10, 0.20, 0.30, 0.40, 0.50],
                   'xticklabels': ['-0.1', '0', '0.1', '0.2', '0.3', '0.4', r'$\geq 0.5$'],
                   'yticks':      [],
                   'xaxstart':    xmin,
                   'xaxend':      xmax,
                   'axoffset':    0,
                   'xlabel':      'Improvement in classification of deleterious alleles (AUROC)',
                   'theme':       'open',
                   'hide':        ['left','right'] }

        if k==0:
            pprops['yticks'] = [0, 10, 20, 30, 40]
            pprops['ylabel'] = 'Number of simulations'
            pprops['hide']   = []
            
        mp.line(             ax=ax_cdel[k], x=[bins+(bin_d/2) for kk in range(len(methods)-1)], y=hist_y, plotprops=lineprops, **pprops)
        mp.plot(type='fill', ax=ax_cdel[k], x=[bins+(bin_d/2) for kk in range(len(methods)-1)], y=hist_y, plotprops=fillprops, **pprops)
            
        #mp.plot(type='bar', ax=ax_cdel[k], x=[bins+(bin_d/2) for kk in range(1, len(methods))], y=hist_y, plotprops=hist_props, **pprops)


    ## c - NRMSE

    for k in range(len(test_sets)):
        xmin, xmax = -0.56, 2.0
        ymin, ymax =  0, 40
        
        dmax = xmax - xmin
        dval = [np.array(y_err[k][i]) - np.array(y_err[k][0]) for i in range(1, len(methods))]
        
        n_bins = N_BINS
        bin_d  = dmax/n_bins
        bins   = np.arange(xmin, xmax+bin_d, bin_d)
        hist_y = []

        for kk in range(len(methods)-1):
            if kk==0:
                hist_y.append([-1000 for i in range(len(bins)+2)])
            else:
                temp_hist_y = [np.sum((bins[i]<=dval[kk]) & (dval[kk]<bins[i+1])) for i in range(len(bins)-1)]
                temp_hist_y.append(np.sum(dval[kk]>=bins[-1]))
                #print(np.sum(temp_hist_y))
                hist_y.append(np.insert(np.append(temp_hist_y, 0), 0, 0))
                
        bins       = np.insert(np.append(bins, xmax+bin_d), 0, 0)
        xmax      += bin_d
        hist_props = dict(lw=AXWIDTH/2, width=0.9*dmax/n_bins, align='center', orientation='vertical', edgecolor=[BKCOLOR], alpha=0.6)
            
        ax_rmse[k].text(xmin/2, 35, 'alternative\nmore\naccurate', color=BKCOLOR, rotation=0, **tprops)
        ax_rmse[k].text(xmax/2, 35, 'MPL more accurate', color=BKCOLOR, rotation=0, **tprops)
        ax_rmse[k].axvline(x=0, **dashlineprops)

        pprops = { 'colors':      colorlist,
                   'xlim':        [xmin, xmax],
                   'ylim':        [ymin, ymax],
                   'xticks':      [-0.5, 0, 0.5, 1.0, 1.5, 2.0],
                   'xticklabels': ['-0.5', '0', '0.5', '1.0', '1.5', r'$\geq 2.0$'],
                   'yticks':      [],
                   'xaxstart':    xmin,
                   'xaxend':      xmax,
                   'axoffset':    0,
                   'xlabel':      'Improvement in error on inferred selection coefficients (NRMSE)',
                   'theme':       'open',
                   'hide':        ['left','right'] }

        if k==0:
            pprops['yticks'] = [0, 10, 20, 30, 40]
            pprops['ylabel'] = 'Number of simulations'
            pprops['hide']   = []
        
        mp.line(             ax=ax_rmse[k], x=[bins+(bin_d/2) for kk in range(len(methods)-1)], y=hist_y, plotprops=lineprops, **pprops)
        mp.plot(type='fill', ax=ax_rmse[k], x=[bins+(bin_d/2) for kk in range(len(methods)-1)], y=hist_y, plotprops=fillprops, **pprops)
        
        #mp.plot(type='bar', ax=ax_rmse[k], x=[bins+(bin_d/2) for kk in range(1, len(methods))], y=hist_y, plotprops=hist_props, **pprops)


    # labels and legend

    labelx = 0.10
    ax_cben[0].text(labelx, box_cben['top'] + 0.01, 'a'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_cdel[0].text(labelx, box_cdel['top'] + 0.01, 'b'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_rmse[0].text(labelx, box_rmse['top'] + 0.01, 'c'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    invt = ax_cben[-1].transData.inverted()
    xy1  = invt.transform((  0, 0))
    xy2  = invt.transform((7.5, 9))
    xy3  = invt.transform((3.0, 9))

    legend_dx1 = xy1[0]-xy2[0]
    legend_dx2 = xy1[0]-xy3[0]
    legend_dy  = xy1[1]-xy2[1]

    legend_x  =  0.6
    legend_d  = -0.015
    legend_y  = 35
    plotprops = DEF_ERRORPROPS.copy()
    plotprops['clip_on'] = False
    for k in range(len(methods)-1):
        mp.error(ax=ax_cben[-1], x=[[legend_x + legend_d]], y=[[legend_y + (k * legend_dy)]],
                 edgecolor=[colorlist[k]], facecolor=[colorbg[k]], plotprops=plotprops, **pprops)
        ax_cben[-1].text(legend_x, legend_y + (k * legend_dy), methods[k+1], ha='left', va='center', **DEF_LABELPROPS)

    # Save figure

    plt.savefig('%s/ed-fig-2-performance%s' % (FIG_DIR, EXT), **FIGPROPS)
    plt.close(fig)

    print('Performance comparison done.')


def plot_supplementary_figure_absolute_delta_s(**pdata):
    """
    Histogram of absolute values of \Delta s for masked variants in each patient/genomic region.
    """
    
    # unpack data
    
    patient_list = pdata['patient_list']
    label_list   = pdata['label_list']
    region_list  = pdata['region_list']
    ds_values    = pdata['ds_values']
    fig_title    = pdata['fig_title']

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 0.94
    goldh   = w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    box_top    = 0.98
    box_bottom = 0.11
    box_left   = 0.10
    box_right  = 0.90
    ddy      = 0.06 / hshrink
    dy       = 0.10 / hshrink

    box_hist = dict(left=box_left, right=box_right, bottom=box_bottom, top=box_top)
    gs_hist  = gridspec.GridSpec(6, 5, wspace=0.15, hspace=0.20, **box_hist)
    i_idx    = 0
    j_idx    = 0

    # iterate through patients/regions and plot |\Delta s| histograms

    for k in range(len(patient_list)):
    
#        ds_max = 1.281
#        #ds_max = 0.006
#        n_bins = 50
        ds_max = 0.40
        n_bins = 41
        bin_ds = ds_max/n_bins
        bins   = np.arange(0, ds_max+bin_ds, bin_ds)
        hist_y = [np.sum((bins[i]<=ds_values[k]) & (ds_values[k]<bins[i+1])) for i in range(len(bins)-1)] #/len(ds_values[k])
        hist_y.append(np.sum(ds_values[k]>=bins[-1]))

        for i in range(len(hist_y)):
            if hist_y[i]>0:
                hist_y[i] = 1 + np.log10(hist_y[i])
            else:
                hist_y[i] = 0

        hist_props = dict(lw=AXWIDTH/2, width=0.9*ds_max/n_bins, align='center', orientation='vertical', edgecolor=[BKCOLOR])

        pprops = { #'xlim':        [0, 1.3],
                   #'xticks':      [ 0, 0.4, 0.8, 1.2, 1.3],
                   #'xticklabels': ['',  '',  '',  '',  ''],
                   'xlim':        [0, 0.41],
                   'xticks':      [ 0, 0.1, 0.2, 0.3, 0.4],
                   'xticklabels': ['',  '',  '',  '',  ''],
                   'ylim':        [0, 4],
                   'yticks':      [],
                   'colors':      [C_MPL],
                   'plotprops':   hist_props,
                   'axoffset':    0.1,
                   'theme':       'open',
                   'hide':        ['left'] }

        if j_idx==0:
            pprops['yticks']      = [0, 1, 2, 3, 4]
            pprops['yticklabels'] = [0, 1, 10, 100, 1000]
            pprops['ylabel']      = 'Counts'
            pprops['hide'].remove('left')
        
        if i_idx==5 or (j_idx>1 and i_idx==4):
            pprops['xticklabels'] = [r'$0$', r'$10$', r'$20$', r'$30$', r'$\geq 40$']
            pprops['xlabel']      = 'Sum of absolute values\nof effects on inferred\nselection coefficients,\n' + r'$\sum_j\|\Delta \hat{s}_{ij}\|$' + ' (%)'

        ax = plt.subplot(gs_hist[i_idx, j_idx])
        mp.plot(type='bar', ax=ax, x=[bins+(bin_ds/2)], y=[hist_y], **pprops)
        
        tprops = dict(ha='center', va='center', family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
        ax.text(0.20, 3.25, label_list[k].upper() + ' ' + (r'$%d\prime$' % int(region_list[k])), **tprops)
        
        if j_idx<4:
            j_idx += 1
        else:
            i_idx += 1
            j_idx  = 0

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, fig_title, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('Delta s histogram done.')


def plot_supplementary_figure_delta_s_correlation(**pdata):
    """
    Scatter plots of \Delta s_ij versus \Delta s_ij for variants in each patient/genomic region.
    """
    
    # unpack data
    
    patient_list = pdata['patient_list']
    region_list  = pdata['region_list']
    ds_values_x  = pdata['ds_values_x']
    ds_values_y  = pdata['ds_values_y']
    fig_title    = pdata['fig_title']

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 0.94
    goldh   = w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    box_top    = 0.98
    box_bottom = 0.11
    box_left   = 0.10
    box_right  = 0.90
    ddy      = 0.06 / hshrink
    dy       = 0.10 / hshrink

    box_hist = dict(left=box_left, right=box_right, bottom=box_bottom, top=box_top)
    gs_hist  = gridspec.GridSpec(6, 5, wspace=0.15, hspace=0.20, **box_hist)
    i_idx    = 0
    j_idx    = 0

    # iterate through patients/regions and plot \Delta s pairs

    for k in range(len(patient_list)):
    
        plotprops = dict(lw=0, marker='o', s=SMALLSIZEDOT, clip_on=False)
        
        pprops = { 'xlim':        [-0.05, 0.05],
                   'xticks':      [-0.05, 0, 0.05],
                   'xticklabels': ['', '', ''],
                   'ylim':        [-0.05, 0.05],
                   'yticks':      [],
                   'yticklabels': [],
                   'colors':      [C_MPL],
                   'plotprops':   plotprops,
                   'axoffset':    0.1,
                   'theme':       'open',
                   'hide':        ['left'] }

        if j_idx==0:
            pprops['yticks']      = [-0.05, 0, 0.05]
            pprops['yticklabels'] = [-5, 0, 5]
            pprops['ylabel']      = r'$\Delta \hat{s}_{ji}$'
            pprops['hide'].remove('left')
        
        if i_idx==5 or (j_idx>1 and i_idx==4):
            pprops['xticklabels'] = [-5, 0, 5]
            pprops['xlabel']      = r'$\Delta \hat{s}_{ij}$' + ', (%)'

        ax = plt.subplot(gs_hist[i_idx, j_idx])
        mp.plot(type='scatter', ax=ax, x=[ds_values_x[k]], y=[ds_values_y[k]], **pprops)
        
        print(st.spearmanr(ds_values_x[k], ds_values_y[k]))
        print(np.min(ds_values_x[k]), np.max(ds_values_x[k]), np.mean(ds_values_x[k]), np.std(ds_values_x[k]))
        print(np.min(ds_values_y[k]), np.max(ds_values_y[k]))
        
        tprops = dict(ha='center', va='center', family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
        ax.text(0.20, 3.25, patient_list[k].lower() + ' ' + (r'$%d\prime$' % int(region_list[k])), **tprops)
        
        if j_idx<4:
            j_idx += 1
        else:
            i_idx += 1
            j_idx  = 0

    # SAVE FIGURE

    plt.savefig('%s/%s.png' % (FIG_DIR, fig_title), dpi = 1200, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('Delta s scatter done.')


def plot_supplementary_figure_delta_s_distance(**pdata):
    """
    Distribution of \Delta s values by genomic distance between the masked and target variants.
    """
    
    # unpack data
    
    ds_values   = np.array(pdata['ds_values'])
    ds_distance = np.array(pdata['ds_distance'])
    fig_title   = pdata['fig_title']

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 0.55
    goldh   = w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    box_top    = 0.95
    box_bottom = 0.12
    box_left   = 0.25
    box_right  = 0.75

    box_dist = dict(left=box_left, right=box_right, bottom=box_bottom, top=box_top)
    gs_dist  = gridspec.GridSpec(1, 1, **box_dist)
    ax_dist  = plt.subplot(gs_dist[0, 0])

    # get \Delta s distribution by distance group
    
#    edges = [0, 100, 700, 1500, 999999]
    edges = [0, 10, 30, 100, 999999]
    group_ds = []
    
    idxs = ds_distance<1
    group_ds.append(ds_values[idxs])
    
    for i in range(len(edges)-1):
        idxs = (ds_distance>edges[i]) & (ds_distance<=edges[i+1])
        group_ds.append(ds_values[idxs])
    
    ds_max = 0.01
    n_bins = 30
    bin_ds = ds_max / n_bins
    bins = np.arange(0, ds_max+bin_ds, bin_ds)
    bin_x = [bins for i in range(len(group_ds))]
    bin_y = []
    y_min = -7

    for i in range(len(group_ds)):
        temp_y = []
        for j in range(len(bins)-1):
            temp_y.append(np.sum((group_ds[i]>=bins[j]) & (group_ds[i]<bins[j+1])))
        temp_y.append(np.sum(group_ds[i]>=bins[-1]))
        temp_y = np.array(temp_y) / np.sum(temp_y)
        for j in range(len(temp_y)):
            if temp_y[j]>0:
                temp_y[j] = np.log10(temp_y[j])
            else:
                temp_y[j] = y_min
        bin_y.append(temp_y)

    # MAKE STEP STYLE HISTOGRAMS FOR EACH DISTANCE GROUP

    lineprops = { 'lw': SIZELINE*2, 'linestyle': '-', 'alpha': 1.0, 'drawstyle': 'steps-mid' }
    #fillprops = { 'lw': 0, 'alpha': 0.2, 'interpolate': True, 'step': 'mid' }

    pprops = { 'xlim':        [-bin_ds, np.max(bins)+bin_ds],
               'xticks':      [0, 0.0025, 0.005, 0.0075, 0.01],
               #'xticklabels': [r'$0$', r'$2.5\times10^{-3}$', r'$5\times10^{-3}$', r'$7.5\times10^{-3}$', r'$\geq 10^{-2}$'],
               'xticklabels': [0, 0.25, 0.5, 0.75, r'$\geq 1$'],
               'ylim':        [-5, 0],
               'yticks':      [-5, -4, -3, -2, -1, 0],
               'yticklabels': ['$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '1'],
               'xlabel':      'Absolute value of effect on inferred selection coefficient, ' + r'$\|\Delta \hat{s}_{ij}\|$' + ' (%)',
               'ylabel':      'Frequency',
               'bins':        bins,
               'combine':     True,
               'plotprops':   lineprops,
               'axoffset':    0.1,
               'theme':       'open' }

    h0 =  0.11
    dh =  0 - 0.11
    l0 =  0.53
    dl =  0.73 - 0.53
    s0 =  1.00
    ds = -1.00

    for k in range(len(group_ds)-1):
        x = k/(len(group_ds) - 1)
        c = hls_to_rgb(h0 + (dh * x), l0 + (dl * x), s0 + (ds * x))
        # Yellow 0.11 0.53 1.00
        # Gray   0.00 0.59 0.00
        mp.line(ax=ax_dist, x=[bin_x[k]], y=[bin_y[k]], colors=[c], zorder=-k, **pprops)

    c = hls_to_rgb(h0 + dh, l0 + dl, s0 + ds)
    mp.plot(type='line', ax=ax_dist, x=[bin_x[-1]], y=[bin_y[-1]], colors=[c], zorder=-len(group_ds), **pprops)

    # Plot legend

    invt = ax_dist.transData.inverted()
    xy1  = invt.transform((  0, 0))
    xy2  = invt.transform((7.5, 9))
    xy3  = invt.transform((3.0, 9))

    legend_dx1 = 1*(xy1[0]-xy2[0])
    legend_dx2 = 1*(xy1[0]-xy3[0])
    legend_dy  = 1*(xy1[1]-xy2[1])

    legend_x =  0.012
    legend_y = -0.5
    legend_d = -0.0005
    legend_t = ['$0$', '$1-10$', '$11-30$', '$31-100$', '$>100$']
    legend_c = [hls_to_rgb(h0 + (dh * x), l0 + (dl * x), s0 + (ds * x)) for x in np.arange(0, 1+1/len(group_ds), 1/len(group_ds))]
    for k in range(len(edges)):
#        legend_t = r'$%d-%d$' % (edges[k], edges[k+1])
#        if k==0:
#            legend_t = r'$%d$' % (edges[k])
#        elif k==len(edges)-2:
#            legend_t = r'$>%d$' % (edges[k])
        mp.line(ax=ax_dist, x=[[legend_x + legend_dx1, legend_x + legend_dx2]],
                y=[[legend_y + (k * legend_dy), legend_y + (k * legend_dy)]],
                colors=[legend_c[k]], plotprops=dict(lw=2*SIZELINE, ls='-', clip_on=False))
        ax_dist.text(legend_x, legend_y + (k * legend_dy), legend_t[k], ha='left', va='center', **DEF_LABELPROPS)
    ax_dist.text(legend_x+0.0005, legend_y - (1.75*legend_dy), 'Distance between variant i\nand target variant j (bp)', ha='center', va='center', **DEF_LABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, fig_title, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('Delta s distance distribution done.')
    

def plot_supplementary_figure_delta_s_icov_distance(**pdata):
    """
    Distribution of \Delta s values by genomic distance between the masked and target variants.
    """
    
    # unpack data
    
    ds_values   = np.array(pdata['ds_values'])
    ds_distance = np.array(pdata['ds_distance'])
    fig_title   = pdata['fig_title']

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 0.55
    goldh   = w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    box_top    = 0.95
    box_bottom = 0.12
    box_left   = 0.08
    box_right  = 0.98
    box_middle = box_left + (box_right - box_left)/2
    box_sep    = 0.06

    box_dist = dict(left=box_left, right=box_middle-box_sep, bottom=box_bottom, top=box_top)
    gs_dist  = gridspec.GridSpec(1, 1, **box_dist)
    ax_dist  = plt.subplot(gs_dist[0, 0])
    
    box_icov = dict(left=box_middle+box_sep, right=box_right, bottom=box_bottom, top=box_top)
    gs_icov  = gridspec.GridSpec(1, 1, **box_icov)
    ax_icov  = plt.subplot(gs_icov[0, 0])
    
    ## a - \Delta s by distance
    
#    edges = [0, 100, 700, 1500, 999999]
    edges = [0, 10, 30, 100, 999999]
    group_ds = []
    
    idxs = ds_distance<1
    group_ds.append(ds_values[idxs])
    
    for i in range(len(edges)-1):
        idxs = (ds_distance>edges[i]) & (ds_distance<=edges[i+1])
        group_ds.append(ds_values[idxs])
    
    ds_max = 0.01
    n_bins = 30
    bin_ds = ds_max / n_bins
    bins = np.arange(0, ds_max+bin_ds, bin_ds)
    bin_x = [bins for i in range(len(group_ds))]
    bin_y = []
    y_min = -7

    for i in range(len(group_ds)):
        temp_y = []
        for j in range(len(bins)-1):
            temp_y.append(np.sum((group_ds[i]>=bins[j]) & (group_ds[i]<bins[j+1])))
        temp_y.append(np.sum(group_ds[i]>=bins[-1]))
        temp_y = np.array(temp_y) / np.sum(temp_y)
        for j in range(len(temp_y)):
            if temp_y[j]>0:
                temp_y[j] = np.log10(temp_y[j])
            else:
                temp_y[j] = y_min
        bin_y.append(temp_y)

    # MAKE STEP STYLE HISTOGRAMS FOR EACH DISTANCE GROUP

    lineprops = { 'lw': SIZELINE*2, 'linestyle': '-', 'alpha': 1.0, 'drawstyle': 'steps-mid' }

    pprops = { 'xlim':        [-bin_ds, np.max(bins)+bin_ds],
               'xticks':      [0, 0.0025, 0.005, 0.0075, 0.01],
               'xticklabels': [0, 0.25, 0.5, 0.75, r'$\geq 1$'],
               'ylim':        [-5, 0],
               'yticks':      [-5, -4, -3, -2, -1, 0],
               'yticklabels': ['$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '1'],
               'xlabel':      'Absolute value of effect on inferred selection coefficient, ' + r'$\|\Delta \hat{s}_{ij}\|$' + ' (%)',
               'ylabel':      'Frequency',
               'bins':        bins,
               'combine':     True,
               'plotprops':   lineprops,
               'axoffset':    0.1,
               'theme':       'open' }

    h0 =  0.11
    dh =  0 - 0.11
    l0 =  0.53
    dl =  0.73 - 0.53
    s0 =  1.00
    ds = -1.00

    for k in range(len(group_ds)-1):
        x = k/(len(group_ds) - 1)
        c = hls_to_rgb(h0 + (dh * x), l0 + (dl * x), s0 + (ds * x))
        # Yellow 0.11 0.53 1.00
        # Gray   0.00 0.59 0.00
        mp.line(ax=ax_dist, x=[bin_x[k]], y=[bin_y[k]], colors=[c], zorder=-k, **pprops)

    c = hls_to_rgb(h0 + dh, l0 + dl, s0 + ds)
    mp.plot(type='line', ax=ax_dist, x=[bin_x[-1]], y=[bin_y[-1]], colors=[c], zorder=-len(group_ds), **pprops)
    
    ax_dist.text(box_left-0.04, box_top+0.02, 'a'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    
    
    ## b - integrated covariance by distance

    edges = [0, 100, 700, 1500, 999999]
    group_icov = []
    
    for i in range(len(edges)-1):
#        group_icov.append(np.fabs(np.loadtxt('%s/binData_%d.txt' % (MLB_DIR, i+1), delimiter=',')))  # Code Ocean
        group_icov.append(np.fabs(np.loadtxt('%s/CatD/binData_%d.txt' % (MLB_DIR, i+1), delimiter=',')))  # GitHub
    
    icov_max = 100
    n_bins   = 30
    bin_icov = icov_max / n_bins
    bins     = np.arange(0, icov_max+bin_icov, bin_icov)
    bin_x    = [bins for i in range(len(group_icov))]
    bin_y    = []
    y_min    = -7

    for i in range(len(group_icov)):
        temp_y = []
        for j in range(len(bins)-1):
            temp_y.append(np.sum((group_icov[i]>=bins[j]) & (group_icov[i]<bins[j+1])))
        temp_y.append(np.sum(group_icov[i]>=bins[-1]))
        temp_y = np.array(temp_y) / np.sum(temp_y)
        for j in range(len(temp_y)):
            if temp_y[j]>0:
                temp_y[j] = np.log10(temp_y[j])
            else:
                temp_y[j] = y_min
        bin_y.append(temp_y)

    lineprops = { 'lw': SIZELINE*2, 'linestyle': '-', 'alpha': 1.0, 'drawstyle': 'steps-mid' }

    pprops = { 'xlim':        [-bin_icov, np.max(bins)+bin_icov],
               'xticks':      [0, 25, 50, 75, 100],
               'xticklabels': [0, 25, 50, 75, r'$\geq 100$'],
               'ylim':        [-5, 0],
               'yticks':      [-5, -4, -3, -2, -1, 0],
               'yticklabels': ['$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '1'],
               'xlabel':      'Absolute value of integrated covariance',
               'ylabel':      'Frequency',
               'bins':        bins,
               'combine':     True,
               'plotprops':   lineprops,
               'axoffset':    0.1,
               'theme':       'open' }

    h0 =  0.11
    dh =  0 - 0.11
    l0 =  0.53
    dl =  0.73 - 0.53
    s0 =  1.00
    ds = -1.00

    for k in range(len(group_icov)-1):
        x = k/(len(group_icov) - 1)
        c = hls_to_rgb(h0 + (dh * x), l0 + (dl * x), s0 + (ds * x))
        # Yellow 0.11 0.53 1.00
        # Gray   0.00 0.59 0.00
        mp.line(ax=ax_icov, x=[bin_x[k]], y=[bin_y[k]], colors=[c], zorder=-k, **pprops)

    c = hls_to_rgb(h0 + dh, l0 + dl, s0 + ds)
    mp.plot(type='line', ax=ax_icov, x=[bin_x[-1]], y=[bin_y[-1]], colors=[c], zorder=-len(group_icov), **pprops)
    
    ax_icov.text(box_middle+box_sep-0.04, box_top+0.02, 'b'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    # Plot legend - a

    invt = ax_dist.transData.inverted()
    xy1  = invt.transform((  0, 0))
    xy2  = invt.transform((7.5, 9))
    xy3  = invt.transform((3.0, 9))

    legend_dx1 = 1*(xy1[0]-xy2[0])
    legend_dx2 = 1*(xy1[0]-xy3[0])
    legend_dy  = 1*(xy1[1]-xy2[1])

    legend_x =  0.0072
    legend_y = -0.5
    legend_d = -0.0005
    legend_t = ['$0$', '$1-10$', '$11-30$', '$31-100$', '$>100$']
    legend_c = [hls_to_rgb(h0 + (dh * x), l0 + (dl * x), s0 + (ds * x)) for x in np.arange(0, 1+1/(len(group_ds)-1), 1/(len(group_ds)-1))] + [c]
    for k in range(len(legend_t)):
        mp.line(ax=ax_dist, x=[[legend_x + legend_dx1, legend_x + legend_dx2]],
                y=[[legend_y + (k * legend_dy), legend_y + (k * legend_dy)]],
                colors=[legend_c[k]], plotprops=dict(lw=2*SIZELINE, ls='-', clip_on=False))
        ax_dist.text(legend_x, legend_y + (k * legend_dy), legend_t[k], ha='left', va='center', **DEF_LABELPROPS)
    ax_dist.text(legend_x+0.0005, legend_y - (1.75*legend_dy), 'Distance between variant i\nand target variant j (bp)', ha='center', va='center', **DEF_LABELPROPS)
    
    # Plot legend - b

    invt = ax_icov.transData.inverted()
    xy1  = invt.transform((  0, 0))
    xy2  = invt.transform((7.5, 9))
    xy3  = invt.transform((3.0, 9))

    legend_dx1 = 1*(xy1[0]-xy2[0])
    legend_dx2 = 1*(xy1[0]-xy3[0])
    legend_dy  = 1*(xy1[1]-xy2[1])

    legend_x =  72
    legend_y = -0.5
    legend_d = -0.0005
    legend_t = ['$1-100$', '$101-700$', '$701-1500$', '$>1500$']
    legend_c = [hls_to_rgb(h0 + (dh * x), l0 + (dl * x), s0 + (ds * x)) for x in np.arange(0, 1+1/(len(group_icov)-1), 1/(len(group_icov)-1))] + [c]
    for k in range(len(legend_t)):
        mp.line(ax=ax_icov, x=[[legend_x + legend_dx1, legend_x + legend_dx2]],
                y=[[legend_y + (k * legend_dy), legend_y + (k * legend_dy)]],
                colors=[legend_c[k]], plotprops=dict(lw=2*SIZELINE, ls='-', clip_on=False))
        ax_icov.text(legend_x, legend_y + (k * legend_dy), legend_t[k], ha='left', va='center', **DEF_LABELPROPS)
    ax_icov.text(legend_x+5, legend_y - (1.75*legend_dy), 'Distance between variant i\nand target variant j (bp)', ha='center', va='center', **DEF_LABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, fig_title, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('Delta s / integrated covariance distance distribution done.')


def plot_figure_delta_s_hive(**pdata):
    """
    Hive plot of \Delta s effects.
    """

    # unpack data
    
    patient_list = pdata['patient_list']
    label_list   = pdata['label_list']
    region_list  = pdata['region_list']
    fig_title    = pdata['fig_title']

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 0.4 # * (0.95 / 0.85)
    goldh   = w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    box_top    = 1.00
    box_bottom = 0.05
    box_left   = 0.025
    box_right  = 0.975

    box_hive = dict(left=box_left, right=box_right, bottom=box_bottom, top=box_top)
    gs_hive  = gridspec.GridSpec(2, 5, wspace=0, hspace=0, **box_hive)
    i_idx    = 0
    j_idx    = 0

    # iterate through patients/regions and plot hive plots

    for k in range(len(patient_list)):
    
        ax = plt.subplot(gs_hive[i_idx, j_idx])
        plot_hive(ax, patient_list[k]+'-'+region_list[k])
        
        tprops = dict(ha='center', va='center', family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
        if i_idx==0:
            ax.text(0, 0.70, label_list[k].upper(), **tprops)
        if j_idx==0:
            ax.text(-1.15, -0.10, r'$%d\prime$' % int(region_list[k]), **tprops)

        if j_idx<2:
            j_idx += 1
        else:
            i_idx += 1
            j_idx  = 0

    # plot legend

    ax_hive = plt.subplot(gs_hive[-2:, -2:])
    
    idx_mask_pos  = [90]
    idx_delta_pos = [90]
    ds_pos        = [ 1]
    idx_mask_neg  = [70]
    idx_delta_neg = [70]
    ds_neg        = [ 1]

    # arc plot for large values of Delta s

    L         = 100
    ds_cutoff = 0.004
    ds_max    = 0.01
    r_min     = 0.025
    r_max     = 0.50
    r_norm    = float(L) / (r_max - r_min)
    arc_mult  = SIZELINE * 3 * ds_max
    arc_alpha = 1

    mask_angle_p =  -np.pi/2
    mask_angle_n = 3*np.pi/2
    ds_pos_angle =   np.pi/6
    ds_neg_angle = 5*np.pi/6
    da           = 0.005

    circ_color  = [hls_to_rgb(0.02, 0.53 * ds + 1. * (1 - ds), 0.83) for ds in ds_pos]
    circ_color += [hls_to_rgb(0.58, 0.53 * ds + 1. * (1 - ds), 0.60) for ds in ds_neg]
    circ_rad    = [[r_min + idx_mask_pos[i]/r_norm, r_min + idx_delta_pos[i]/r_norm] for i in range(len(ds_pos))]
    circ_rad   += [[r_min + idx_mask_neg[i]/r_norm, r_min + idx_delta_neg[i]/r_norm] for i in range(len(ds_neg))]
    circ_arc    = [dict(lw=arc_mult * np.fabs(0.1) / ds_cutoff, alpha=arc_alpha) for ds in ds_pos]
    circ_arc   += [dict(lw=arc_mult * np.fabs(0.1) / ds_cutoff, alpha=arc_alpha) for ds in ds_neg]
    circ_angle  = [[mask_angle_p+da/r[0], ds_pos_angle-da/r[1]] for r in circ_rad[:len(ds_pos)]]
    circ_angle += [[mask_angle_n-da/r[0], ds_neg_angle+da/r[1]] for r in circ_rad[len(ds_pos):]]
    circ_bezrad = [(r[0]+r[1])/2 for r in circ_rad]
    circ_x      = [i for i in ds_pos] + [i for i in ds_neg]
    circ_y      = [i for i in ds_pos] + [i for i in ds_neg]
    
    # CD8+ T cell epitope legend

    r_mid = r_min + (20 / r_norm)
    scatter_x = [r_mid * np.cos(mask_angle_p),
                 r_mid * np.cos(ds_pos_angle),
                 r_mid * np.cos(ds_neg_angle) ]
    scatter_y = [r_mid * np.sin(mask_angle_p),
                 r_mid * np.sin(ds_pos_angle),
                 r_mid * np.sin(ds_neg_angle) ]
    smallprops = dict(lw=0, marker='o', s=1.0*SMALLSIZEDOT, zorder=9999, clip_on=False)
    bigprops   = dict(lw=0, marker='o', s=1.7*SMALLSIZEDOT, zorder=9998, clip_on=False)
    mp.scatter(ax=ax_hive, x=[scatter_x], y=[scatter_y], colors=['#FFFFFF'], plotprops=smallprops)
    mp.scatter(ax=ax_hive, x=[scatter_x], y=[scatter_y], colors=[BKCOLOR],   plotprops=  bigprops)
    
    # add legend labels

    ## epitope label
    label_r  = r_mid
    ddx, ddy = 0.09/2, 0.15/2
    label_x  = [label_r * np.cos(ds_neg_angle)]
    label_y  = [label_r * np.sin(ds_neg_angle) + 0.5*ddy]
    label_x  = label_x + [label_x[0]]
    label_y  = label_y + [label_y[0] + 2*ddy]

    txtprops  = dict(ha='center', va='bottom', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
    plotprops = dict(lw=AXWIDTH, ls='-', clip_on=False)
    
    ax_hive.text(label_x[-1], label_y[-1] + 0.125*ddy, 'CD8+ T cell\nepitope', **txtprops)
    mp.line(ax=ax_hive, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)
    
    ## 5' -- 3' label
    arrow_x  = np.array([r_min * np.cos(ds_pos_angle), r_max * np.cos(ds_pos_angle)]) - 0.03/2
    arrow_y  = np.array([r_min * np.sin(ds_pos_angle), r_max * np.sin(ds_pos_angle)]) + 0.10/2
    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=180*ds_pos_angle/np.pi)
    
    ddr = 0.05/2
    ax_hive.text(arrow_x[0] + 1.0 * ddr * np.cos(ds_pos_angle), arrow_y[0] + 0.1/2 + 1.0 * ddr * np.sin(ds_pos_angle), r'$5\prime$', **txtprops)
    ax_hive.text(arrow_x[1] - 2.5 * ddr * np.cos(ds_pos_angle), arrow_y[1] + 0.1/2 - 2.5 * ddr * np.sin(ds_pos_angle), r'$3\prime$', **txtprops)

    ax_hive.annotate('', xy=(arrow_x[1], arrow_y[1]), xytext=(arrow_x[0], arrow_y[0]),
                     arrowprops=dict(color=BKCOLOR, lw=0, shrink=0.0, width=AXWIDTH, headwidth=AXWIDTH*5, headlength=AXWIDTH*5),)
    
    ## focus/target labels
    label_r = 1.15/2
    tprops  = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_hive.text(label_r * np.cos(mask_angle_n), label_r * np.sin(mask_angle_n), 'Variant i', **tprops)
    label_r = 1.1/2
    tprops['ha'] = 'left'
    ax_hive.text(label_r * np.cos(ds_pos_angle), label_r * np.sin(ds_pos_angle), 'Target\nvariant j', **tprops)
    tprops['ha'] = 'right'
    ax_hive.text(label_r * np.cos(ds_neg_angle), label_r * np.sin(ds_neg_angle), 'Target\nvariant j', **tprops)
    
    ## negative link label
    label_r  = r_min + idx_mask_neg[0]/r_norm - 0.11/2
    ddx, ddy = 0.09/2, 0.15/2
    label_x  = [label_r * np.cos(np.pi)]
    label_y  = [label_r * np.sin(np.pi) - 0.5*ddy]
    label_x  = label_x + [label_x[0] - ddx, label_x[0] - 1.5*ddx]
    label_y  = label_y + [label_y[0] - ddy, label_y[0] -     ddy]
    
    txtprops  = dict(ha='right', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
    plotprops = dict(lw=AXWIDTH, ls='-', clip_on=False)
    
    ax_hive.text(label_x[-1] - 0.5*ddx, label_y[-1], 'Variant i decreases\ntarget ' + r'$\hat{s}_j$' + ', ' + r'$\Delta \hat{s}_{ij}<0$', **txtprops)
    mp.line(ax=ax_hive, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)
    
    ## positive link label
    label_r  = r_min + idx_mask_pos[0]/r_norm - 0.15/2
    ddx, ddy = 0.09/2, 0.15/2
    label_x  = [label_r * np.cos(0)]
    label_y  = [label_r * np.sin(0) - 0.5*ddy]
    label_x  = label_x + [label_x[0] + ddx, label_x[0] + 1.5*ddx]
    label_y  = label_y + [label_y[0] - ddy, label_y[0] -     ddy]
    
    txtprops  = dict(ha='left', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
    plotprops = dict(lw=AXWIDTH, ls='-', clip_on=False)
    
    ax_hive.text(label_x[-1] + 0.5*ddx, label_y[-1], 'Variant i increases\ntarget ' + r'$\hat{s}_j$' + ', ' + r'$\Delta \hat{s}_{ij}>0$', **txtprops)
    mp.line(ax=ax_hive, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)
    
    ## color bar
    ddx      =  0 #0.525
    ddy      = -0.15
    label_x  =  0
    label_y  = -0.96
    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_hive.text(label_x+ddx, label_y+ddy, 'Effect on inferred selection\ncoefficient, ' + r'$\Delta \hat{s}_{ij}$' + ' (%)', **txtprops)
    
    dbox      = 0.05
    rec_props = dict(height=dbox, width=dbox, ec=None, lw=AXWIDTH/2, clip_on=False)
    for i in range(-5, 5+1, 1):
        c = BKCOLOR
        t = i/5
        if t>0:
            c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
        else:
            c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
        rec = matplotlib.patches.Rectangle(xy=((i-0.5)*dbox + ddx, -0.75 + ddy), fc=c, **rec_props)
        ax_hive.add_artist(rec)

    txtprops = dict(ha='center', va='top', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_hive.text(-5.25*dbox + ddx, -0.78 + ddy, -1, **txtprops)
    ax_hive.text(             ddx, -0.78 + ddy,  0, **txtprops)
    ax_hive.text( 5.00*dbox + ddx, -0.78 + ddy,  1, **txtprops)
    
    # Make plot

    pprops = { 'colors':   circ_color,
               'xlim':     [-1.0, 1.0],
               'ylim':     [-1.3, 0.7],
               'size':     L,
               'rad':      circ_rad,
               'arcprops': circ_arc,
               'angle':    circ_angle,
               'bezrad':   circ_bezrad,
               'noaxes':   True }

    plotprops = dict(lw=AXWIDTH, ls='-', zorder=999)
    line_x = [[r_min * np.cos(mask_angle_p), r_max * np.cos(mask_angle_p)],
              [r_min * np.cos(ds_pos_angle), r_max * np.cos(ds_pos_angle)],
              [r_min * np.cos(ds_neg_angle), r_max * np.cos(ds_neg_angle)] ]
    line_y = [[r_min * np.sin(mask_angle_p), r_max * np.sin(mask_angle_p)],
              [r_min * np.sin(ds_pos_angle), r_max * np.sin(ds_pos_angle)],
              [r_min * np.sin(ds_neg_angle), r_max * np.sin(ds_neg_angle)] ]
    mp.line(ax=ax_hive, x=line_x, y=line_y, colors=[BKCOLOR for i in range(len(line_x))], plotprops=plotprops)

    mp.plot(type='circos', ax=ax_hive, x=circ_x, y=circ_y, **pprops)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, fig_title, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('%s done.' % fig_title)


def plot_supplementary_figure_delta_s_hive(**pdata):
    """
    Hive plot of \Delta s effects.
    """

    # unpack data
    
    patient_list = pdata['patient_list']
    label_list   = pdata['label_list']
    region_list  = pdata['region_list']
    fig_title    = pdata['fig_title']

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 1.2
    goldh   = w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    box_top    = 1.00
    box_bottom = 0.05
    box_left   = 0.025
    box_right  = 0.975

    box_hive = dict(left=box_left, right=box_right, bottom=box_bottom, top=box_top)
    gs_hive  = gridspec.GridSpec(6, 5, wspace=0, hspace=0, **box_hive)
    i_idx    = 0
    j_idx    = 0

    # iterate through patients/regions and plot hive plots

    for k in range(len(patient_list)):
    
        ax = plt.subplot(gs_hive[i_idx, j_idx])
        plot_hive(ax, patient_list[k]+'-'+region_list[k])
        
        tprops = dict(ha='center', va='center', family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
        ax.text(0, -1.15, label_list[k].upper() + ' ' + (r'$%d\prime$' % int(region_list[k])), **tprops)

        if i_idx<4 and j_idx<4:
            j_idx += 1
        elif i_idx>=4 and j_idx<2:
            j_idx += 1
        else:
            i_idx += 1
            j_idx  = 0

    # plot legend

    ax_hive = plt.subplot(gs_hive[-2:, -2:])
    
    idx_mask_pos  = [90]
    idx_delta_pos = [90]
    ds_pos        = [ 1]
    idx_mask_neg  = [70]
    idx_delta_neg = [70]
    ds_neg        = [ 1]

    # arc plot for large values of Delta s

    L         = 100
    ds_cutoff = 0.004
    ds_max    = 0.01
    r_min     = 0.025
    r_max     = 0.50
    r_norm    = float(L) / (r_max - r_min)
    arc_mult  = SIZELINE * 3 * ds_max
    arc_alpha = 1

    mask_angle_p =  -np.pi/2
    mask_angle_n = 3*np.pi/2
    ds_pos_angle =   np.pi/6
    ds_neg_angle = 5*np.pi/6
    da           = 0.005

    circ_color  = [hls_to_rgb(0.02, 0.53 * ds + 1. * (1 - ds), 0.83) for ds in ds_pos]
    circ_color += [hls_to_rgb(0.58, 0.53 * ds + 1. * (1 - ds), 0.60) for ds in ds_neg]
    circ_rad    = [[r_min + idx_mask_pos[i]/r_norm, r_min + idx_delta_pos[i]/r_norm] for i in range(len(ds_pos))]
    circ_rad   += [[r_min + idx_mask_neg[i]/r_norm, r_min + idx_delta_neg[i]/r_norm] for i in range(len(ds_neg))]
    circ_arc    = [dict(lw=arc_mult * np.fabs(0.1) / ds_cutoff, alpha=arc_alpha) for ds in ds_pos]
    circ_arc   += [dict(lw=arc_mult * np.fabs(0.1) / ds_cutoff, alpha=arc_alpha) for ds in ds_neg]
    circ_angle  = [[mask_angle_p+da/r[0], ds_pos_angle-da/r[1]] for r in circ_rad[:len(ds_pos)]]
    circ_angle += [[mask_angle_n-da/r[0], ds_neg_angle+da/r[1]] for r in circ_rad[len(ds_pos):]]
    circ_bezrad = [(r[0]+r[1])/2 for r in circ_rad]
    circ_x      = [i for i in ds_pos] + [i for i in ds_neg]
    circ_y      = [i for i in ds_pos] + [i for i in ds_neg]
    
    # CD8+ T cell epitope legend

    r_mid = r_min + (20 / r_norm)
    scatter_x = [r_mid * np.cos(mask_angle_p),
                 r_mid * np.cos(ds_pos_angle),
                 r_mid * np.cos(ds_neg_angle) ]
    scatter_y = [r_mid * np.sin(mask_angle_p),
                 r_mid * np.sin(ds_pos_angle),
                 r_mid * np.sin(ds_neg_angle) ]
    smallprops = dict(lw=0, marker='o', s=1.0*SMALLSIZEDOT, zorder=9999, clip_on=False)
    bigprops   = dict(lw=0, marker='o', s=1.7*SMALLSIZEDOT, zorder=9998, clip_on=False)
    mp.scatter(ax=ax_hive, x=[scatter_x], y=[scatter_y], colors=['#FFFFFF'], plotprops=smallprops)
    mp.scatter(ax=ax_hive, x=[scatter_x], y=[scatter_y], colors=[BKCOLOR],   plotprops=  bigprops)
    
    # add legend labels

    ## epitope label
    label_r  = r_mid
    ddx, ddy = 0.09/2, 0.15/2
    label_x  = [label_r * np.cos(ds_neg_angle)]
    label_y  = [label_r * np.sin(ds_neg_angle) + 0.5*ddy]
    label_x  = label_x + [label_x[0]]
    label_y  = label_y + [label_y[0] + 2*ddy]

    txtprops  = dict(ha='center', va='bottom', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
    plotprops = dict(lw=AXWIDTH, ls='-', clip_on=False)
    
    ax_hive.text(label_x[-1], label_y[-1] + 0.125*ddy, 'CD8+ T cell\nepitope', **txtprops)
    mp.line(ax=ax_hive, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)
    
    ## 5' -- 3' label
    arrow_x  = np.array([r_min * np.cos(ds_pos_angle), r_max * np.cos(ds_pos_angle)]) - 0.03/2
    arrow_y  = np.array([r_min * np.sin(ds_pos_angle), r_max * np.sin(ds_pos_angle)]) + 0.10/2
    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=180*ds_pos_angle/np.pi)
    
    ddr = 0.05/2
    ax_hive.text(arrow_x[0] + 1.0 * ddr * np.cos(ds_pos_angle), arrow_y[0] + 0.1/2 + 1.0 * ddr * np.sin(ds_pos_angle), r'$5\prime$', **txtprops)
    ax_hive.text(arrow_x[1] - 2.5 * ddr * np.cos(ds_pos_angle), arrow_y[1] + 0.1/2 - 2.5 * ddr * np.sin(ds_pos_angle), r'$3\prime$', **txtprops)

    ax_hive.annotate('', xy=(arrow_x[1], arrow_y[1]), xytext=(arrow_x[0], arrow_y[0]),
                     arrowprops=dict(color=BKCOLOR, lw=0, shrink=0.0, width=AXWIDTH, headwidth=AXWIDTH*5, headlength=AXWIDTH*5),)
    
    ## focus/target labels
    label_r = 1.15/2
    tprops  = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_hive.text(label_r * np.cos(mask_angle_n), label_r * np.sin(mask_angle_n), 'Variant i', **tprops)
    label_r = 1.1/2
    tprops['ha'] = 'left'
    ax_hive.text(label_r * np.cos(ds_pos_angle), label_r * np.sin(ds_pos_angle), 'Target\nvariant j', **tprops)
    tprops['ha'] = 'right'
    ax_hive.text(label_r * np.cos(ds_neg_angle), label_r * np.sin(ds_neg_angle), 'Target\nvariant j', **tprops)
    
    ## negative link label
    label_r  = r_min + idx_mask_neg[0]/r_norm - 0.11/2
    ddx, ddy = 0.09/2, 0.15/2
    label_x  = [label_r * np.cos(np.pi)]
    label_y  = [label_r * np.sin(np.pi) - 0.5*ddy]
    label_x  = label_x + [label_x[0] - ddx, label_x[0] - 1.5*ddx]
    label_y  = label_y + [label_y[0] - ddy, label_y[0] -     ddy]
    
    txtprops  = dict(ha='right', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
    plotprops = dict(lw=AXWIDTH, ls='-', clip_on=False)
    
    ax_hive.text(label_x[-1] - 0.5*ddx, label_y[-1], 'Variant i decreases\ntarget ' + r'$\hat{s}_j$' + ', ' + r'$\Delta \hat{s}_{ij}<0$', **txtprops)
    mp.line(ax=ax_hive, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)
    
    ## positive link label
    label_r  = r_min + idx_mask_pos[0]/r_norm - 0.15/2
    ddx, ddy = 0.09/2, 0.15/2
    label_x  = [label_r * np.cos(0)]
    label_y  = [label_r * np.sin(0) - 0.5*ddy]
    label_x  = label_x + [label_x[0] + ddx, label_x[0] + 1.5*ddx]
    label_y  = label_y + [label_y[0] - ddy, label_y[0] -     ddy]
    
    txtprops  = dict(ha='left', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
    plotprops = dict(lw=AXWIDTH, ls='-', clip_on=False)
    
    ax_hive.text(label_x[-1] + 0.5*ddx, label_y[-1], 'Variant i increases\ntarget ' + r'$\hat{s}_j$' + ', ' + r'$\Delta \hat{s}_{ij}>0$', **txtprops)
    mp.line(ax=ax_hive, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)
    
    ## color bar
    ddx      =  0.425
    ddy      = -0.10
    label_x  =  0
    label_y  = -0.96
    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_hive.text(label_x+ddx, label_y+ddy, 'Effect on inferred selection\ncoefficient, ' + r'$\Delta \hat{s}_{ij}$' + ' (%)', **txtprops)
    
    dbox      = 0.05
    rec_props = dict(height=dbox, width=dbox, ec=None, lw=AXWIDTH/2, clip_on=False)
    for i in range(-5, 5+1, 1):
        c = BKCOLOR
        t = i/5
        if t>0:
            c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
        else:
            c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
        rec = matplotlib.patches.Rectangle(xy=((i-0.5)*dbox + ddx, -0.75 + ddy), fc=c, **rec_props)
        ax_hive.add_artist(rec)

    txtprops = dict(ha='center', va='top', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_hive.text(-5.25*dbox + ddx, -0.78 + ddy, -1, **txtprops)
    ax_hive.text(             ddx, -0.78 + ddy,  0, **txtprops)
    ax_hive.text( 5.00*dbox + ddx, -0.78 + ddy,  1, **txtprops)
    
    # Make plot

    pprops = { 'colors':   circ_color,
               'xlim':     [-1.0, 1.0],
               'ylim':     [-1.0, 1.0],
               'size':     L,
               'rad':      circ_rad,
               'arcprops': circ_arc,
               'angle':    circ_angle,
               'bezrad':   circ_bezrad,
               'noaxes':   True }

    plotprops = dict(lw=AXWIDTH, ls='-', zorder=999)
    line_x = [[r_min * np.cos(mask_angle_p), r_max * np.cos(mask_angle_p)],
              [r_min * np.cos(ds_pos_angle), r_max * np.cos(ds_pos_angle)],
              [r_min * np.cos(ds_neg_angle), r_max * np.cos(ds_neg_angle)] ]
    line_y = [[r_min * np.sin(mask_angle_p), r_max * np.sin(mask_angle_p)],
              [r_min * np.sin(ds_pos_angle), r_max * np.sin(ds_pos_angle)],
              [r_min * np.sin(ds_neg_angle), r_max * np.sin(ds_neg_angle)] ]
    mp.line(ax=ax_hive, x=line_x, y=line_y, colors=[BKCOLOR for i in range(len(line_x))], plotprops=plotprops)

    mp.plot(type='circos', ax=ax_hive, x=circ_x, y=circ_y, **pprops)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, fig_title, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('Hive plots done.')


def plot_supplementary_figure_max_dx(**pdata):
    """
    Distribution of maximum observed change in variant frequency per generation.
    """
    
    # unpack data
    
    dx_sweep  = np.array(pdata['dx_sweep'])
    dx_other  = np.array(pdata['dx_other'])
    fig_title = pdata['fig_title']

    # PLOT FIGURE
    
    # SHOW FREQUENCY INSTEAD

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 0.55
    goldh   = w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    box_top    = 0.95
    box_bottom = 0.12
    box_left   = 0.25
    box_right  = 0.75

    box_dist = dict(left=box_left, right=box_right, bottom=box_bottom, top=box_top)
    gs_dist  = gridspec.GridSpec(1, 1, **box_dist)
    ax_dist  = plt.subplot(gs_dist[0, 0])

    # get dx bins
    
    dx_max  = 0.143
    n_bins  = 20
    bin_dx  = dx_max/n_bins
    bins    = np.arange(0, dx_max+bin_dx, bin_dx)
    hist_ys = []
    for dx_values in [dx_sweep, dx_other]:
        dx_values = np.array(dx_values)
        hist_y = [np.sum((bins[i]<=dx_values) & (dx_values<bins[i+1])) for i in range(len(bins)-1)]
        hist_y.append(np.sum(dx_values>=bins[-1]))
        hist_ys.append(np.array(hist_y)/len(dx_values))

    for k in range(len(hist_ys)):
        for i in range(len(hist_ys[k])):
            if hist_ys[k][i]>0:
                hist_ys[k][i] = np.log10(hist_ys[k][i])
            else:
                hist_ys[k][i] = -10

    # make figure

    lineprops = { 'lw': SIZELINE*2, 'linestyle': '-', 'alpha': 1.0, 'drawstyle': 'steps-mid' }
    hist_props = dict(lw=AXWIDTH/2, width=0.9*dx_max/n_bins, align='center',
                      orientation='vertical', edgecolor=[BKCOLOR])

    pprops = { 'xlim':        [-bin_dx, np.max(bins)+bin_dx],
               'xticks':      [ 0, 0.03, 0.06, 0.09, 0.12, 0.15],
               'xticklabels': [ 0,    3,    6,    9,   12,   15],
               'ylim':        [-4,  0],
               'yticks':      [-4, -3, -2, -1, 0],
               'yticklabels': ['$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '1'],
               'xlabel':      'Maximum observed change in variant frequency per day (%)',
               'ylabel':      'Frequency',
               #'colors':      [C_MPL],
               'plotprops':   lineprops,
               'axoffset':    0.1,
               'theme':       'open' }

#    mp.plot(type='bar', ax=ax, x=[bins+(bin_ds/2), bins+(bin_ds/2)], y=hist_ys, **pprops)

    h0 =  0.11
    dh =  0 - 0.11
    l0 =  0.53
    dl =  0.73 - 0.53
    s0 =  1.00
    ds = -1.00

    c = hls_to_rgb(h0, l0, s0)
    mp.line(ax=ax_dist, x=[bins], y=[hist_ys[0]], colors=[c], **pprops)

    c = hls_to_rgb(h0 + dh, l0 + dl, s0 + ds)
    mp.plot(type='line', ax=ax_dist, x=[bins], y=[hist_ys[-1]], colors=[c], **pprops)

    # Plot legend

    invt = ax_dist.transData.inverted()
    xy1  = invt.transform((  0, 0))
    xy2  = invt.transform((7.5, 9))
    xy3  = invt.transform((3.0, 9))

    legend_dx1 = 1*(xy1[0]-xy2[0])
    legend_dx2 = 1*(xy1[0]-xy3[0])
    legend_dy  = 1.2*(xy1[1]-xy2[1])

    legend_x =  0.17
    legend_y = -0.5
    legend_d = -0.0005
    legend_t = ['Highly influential variants', 'Other variants']
    legend_c = [hls_to_rgb(h0, l0, s0), hls_to_rgb(h0 + dh, l0 + dl, s0 + ds)]
    for k in range(len(legend_t)):
        mp.line(ax=ax_dist, x=[[legend_x + legend_dx1, legend_x + legend_dx2]],
                y=[[legend_y + (1.0 * k * legend_dy), legend_y + (1.0 * k * legend_dy)]],
                colors=[legend_c[k]], plotprops=dict(lw=2*SIZELINE, ls='-', clip_on=False))
        ax_dist.text(legend_x, legend_y + (1.0 * k * legend_dy), legend_t[k], ha='left', va='center', **DEF_LABELPROPS)
    #ax_dist.text(legend_x+0.0005, legend_y - (1.75*legend_dy), 'Distance between variant i\nand target variant j (bp)', ha='center', va='center', **DEF_LABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, fig_title, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('Delta x distribution done.')


def plot_supplementary_figure_epitope(**pdata):
    """
    Alternate epitope escape mutation frequencies, inferred selection coefficients, and linkage.
    """

    # unpack data
    
    patient       = pdata['patient']
    region        = pdata['region']
    epitope       = pdata['epitope']
    epitope_range = pdata['epitope_range']
    epitope_label = pdata['epitope_label']
    cov_label     = pdata['cov_label']
    label2ddr     = pdata['label2ddr']
    legend_loc    = pdata['legend_loc']
    traj_ticks    = pdata['traj_ticks']
    sel_ticks     = pdata['sel_ticks']
    sel_minors    = pdata['sel_minors']
    sel_space     = pdata['sel_space']
    fig_title     = pdata['fig_title']
    tag           = patient+'-'+region
    
    # process stored data
    
    df_poly = pd.read_csv('%s/analysis/%s-poly.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    df_esc  = df_poly[(df_poly.epitope==epitope) & (df_poly.nucleotide!=df_poly.TF)]

    times = [int(i.split('_')[-1]) for i in df_esc.columns if 'f_at_' in i]
    times.sort()
    
    var_tag  = []
    var_smpl = []
    var_sind = []
    var_traj = []
    curr_HXB2 = 8988
    #curr_aln  = 4015
    #curr_char = 'a'
    for df_iter, df_entry in df_esc.iterrows():
        if df_entry.nucleotide=='-':
            continue
        var_tag.append(str(df_entry.HXB2_index)+df_entry.nucleotide)
#        if pd.notnull(df_entry.HXB2_index):
#            var_tag.append(str(int(df_entry.HXB2_index))+df_entry.nucleotide)
#            curr_HXB2 = int(df_entry.HXB2_index)
#            curr_aln  = int(df_entry.alignment_index)
#        else:
#            if int(df_entry.alignment_index)!=curr_aln:
#                curr_aln  = int(df_entry.alignment_index)
#                curr_char = chr(ord(curr_char) + 1)
#            var_tag.append(str(int(curr_HXB2))+curr_char+df_entry.nucleotide)
        var_traj.append([df_entry['f_at_%d' % t] for t in times])
        var_smpl.append(df_entry.s_MPL)
        var_sind.append(df_entry.s_SL)

    var_c = sns.husl_palette(len(var_traj))

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 0.93
    goldh   = 1 * w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    box_top  = 0.95
    dy       = 0.10
    y0       = 0.09/hshrink
    y1       = 0.08/hshrink
    y2       = 0.50/hshrink
    s_spc    = np.max([5, len(var_c)]) * 0.02
    box_traj = dict(left=0.30,                 right=0.70,                     bottom=box_top-y0,                top=box_top)
    box_smpl = dict(left=0.30,                 right=0.30+s_spc,               bottom=box_top-y0-y1-(0.8*dy),    top=box_top-y0-(0.8*dy))
    box_sind = dict(left=0.30+s_spc+sel_space, right=0.30+(2*s_spc)+sel_space, bottom=box_top-y0-y1-(0.8*dy),    top=box_top-y0-(0.8*dy))
    box_circ = dict(left=0.25,                 right=0.75,                     bottom=box_top-y0-y1-y2-(1.9*dy), top=box_top-y0-y1-(1.9*dy))

    gs_traj = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_traj)
    gs_smpl = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_smpl)
    gs_sind = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_sind)
    gs_circ = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_circ)
    ax_traj = plt.subplot(gs_traj[0, 0])
    ax_smpl = plt.subplot(gs_smpl[0, 0])
    ax_sind = plt.subplot(gs_sind[0, 0])
    ax_circ = plt.subplot(gs_circ[0, 0])

    dy = 0.02/hshrink

    ## a -- trajectory plot

    pprops = { 'xticks':      traj_ticks,
               'yticks':      [0, 1],
               'yminorticks': [0.25, 0.5, 0.75],
               'nudgey':      1.1,
               'xlabel':      'Time (days)',
               'ylabel':      'Variant frequency\nin %s epitope\n' % cov_label,
               'plotprops':   {'lw': SIZELINE*1.5, 'ls': '-', 'alpha': 1.0 },
               'axoffset':    0.1,
               'theme':       'open' }

    for i in range(len(var_tag)-1):
        xdat = [times]
        ydat = [var_traj[i]]
        mp.line(ax=ax_traj, x=xdat, y=ydat, colors=[var_c[i]], **pprops)

    xdat = [times]
    ydat = [var_traj[len(var_tag)-1]]
    mp.plot(type='line', ax=ax_traj, x=xdat, y=ydat, colors=[var_c[len(var_tag)-1]], **pprops)

    ax_traj.text(0.25, box_traj['top']+dy, 'a'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## b1 -- selection coefficients inferred by MPL

    xmin, xmax = 0, len(var_tag)
    diff = 5. - len(var_c)
    if diff>0:
        xmin -= diff/2
        xmax += diff/2

    hist_props = dict(lw=SIZELINE/2, width=0.6, align='center', orientation='vertical',
                      edgecolor=[BKCOLOR for i in range(len(var_tag))])

    bar_x  = [i+0.5 for i in range(len(var_tag))]
    pprops = { 'colors':      [var_c],
               'xlim':        [xmin, xmax],
               'ylim':        [np.min(sel_ticks), np.max(sel_ticks)],
               'xticks':      bar_x,
               'xticklabels': var_tag,
               'yticks':      sel_ticks,
               'yticklabels': [int(100*l) for l in sel_ticks],
               'yminorticks': sel_minors,
               'ylabel':      'Inferred selection\ncoefficient, '+r'$\hat{s}$'+' (%)',
               'theme':       'open',
               'hide' :       [] }

    mp.plot(type='bar', ax=ax_smpl, x=[bar_x], y=[var_smpl], plotprops=hist_props, **pprops)
    plt.setp(ax_smpl.xaxis.get_majorticklabels(), rotation=90)

    transFigureInv = fig.transFigure.inverted()
    labelprops     = dict(color=BKCOLOR, ha='center', va='top', family=FONTFAMILY, size=SIZELABEL,
                          clip_on=False, transform=fig.transFigure)
    ax_smpl.text(box_smpl['left']+(box_smpl['right']-box_smpl['left'])/2, box_smpl['top'], 'MPL', **labelprops)

    ax_smpl.text(0.25, box_smpl['top']+dy, 'b'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## b2 -- selection coefficients inferred by SL

    pprops = { 'colors':      [var_c],
               'xlim':        [xmin, xmax],
               'ylim':        [np.min(sel_ticks), np.max(sel_ticks)],
               'xticks':      bar_x,
               'xticklabels': var_tag,
               'yticks':      [],
               'theme':       'open',
               'hide' :       ['left','right'] }

    mp.plot(type='bar', ax=ax_sind, x=[bar_x], y=[var_sind], plotprops=hist_props, **pprops)
    plt.setp(ax_sind.xaxis.get_majorticklabels(), rotation=90)

    ax_sind.text(box_sind['left']+(box_sind['right']-box_sind['left'])/2, box_sind['top'], 'Independent\nmodel', **labelprops)

    # add background

    cBG = '#F5F5F5'
    bg  = ax_sind.axis()
    ddx = 0.01
    ddy = 0.01
    rec = matplotlib.patches.Rectangle(xy=(box_sind['left']-(0.1*ddx), box_sind['bottom']-(0.7*ddy)), width=box_sind['right']-box_sind['left']+(0.2*ddx),
                                       height=box_sind['top']-box_sind['bottom']+(1.7*ddy), transform=fig.transFigure, ec=None, fc=cBG, clip_on=False, zorder=-100)
    rec = ax_sind.add_patch(rec)

    ## c -- circle plot

    sig_s, sig_site_real, sig_nuc_idx, epitope_start, epitope_end = plot_circle(ax_circ, tag, epitope_range, epitope_label, cov_label, label2ddr)
    
    ax_circ.text(0.25, box_circ['top'], 'c'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    # MAKE LEGEND

    invt = ax_circ.transData.inverted()
    xy1  = invt.transform((0,0))
    xy2  = invt.transform((0,9))
    
    coef_legend_x  = 1.30
    coef_legend_dx = 0.05
    coef_legend_y  = 1.00
    coef_legend_dy = xy1[1]-xy2[1]

    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    s_mult   = 20*SMALLSIZEDOT
    ex_s     = [-0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.10]
    show_s   = [   1,     0,     0,     0,     0, 1,    0,    0,    0,    0,    1]
    c_s      = [C_DEL, C_BEN]
    for i in range(len(ex_s)):
        plotprops = dict(lw=0, marker='o', s=np.fabs(ex_s[i])*s_mult, clip_on=False)
        mp.scatter(ax=ax_circ, x=[[coef_legend_x + i*coef_legend_dx]], y=[[coef_legend_y]], colors=[c_s[ex_s[i]>0]], plotprops=plotprops)
        if show_s[i]:
            ax_circ.text(coef_legend_x + i*coef_legend_dx, coef_legend_y + 0.75*coef_legend_dy, '%d' % (100*ex_s[i]), **txtprops)
    ax_circ.text(coef_legend_x + 5*coef_legend_dx, coef_legend_y + (2.25*coef_legend_dy), 'Inferred selection\ncoefficient, $\hat{s}$ (%)', **txtprops)

    coef_legend_x = 1.30
    coef_legend_y = coef_legend_y + 5.5 * coef_legend_dy
    txtprops = dict(ha='left', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    arc_props = dict(lw=AXWIDTH/2, ec=BKCOLOR, fc='#f2f2f2', alpha=0.8, clip_on=False)
    epatch    = matplotlib.patches.Rectangle(xy=(coef_legend_x-0.025, coef_legend_y+(0.0*coef_legend_dy)), width=0.05, height=-1.*coef_legend_dy, **arc_props)
    ax_circ.add_artist(epatch)
    ax_circ.text(coef_legend_x+(1*coef_legend_dx), coef_legend_y-(coef_legend_dy/2), 'CD8+ T cell\nepitope', **txtprops)

    # SAVE FIGURE

    plt.savefig('%s/%s%s' % (FIG_DIR, fig_title, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('%s done.' % fig_title)


def plot_supplementary_figure_cap256_recombination(**pdata):
    """
    Plot patterns of recombination in CAP256.
    """

    # unpack data
    
    patient = pdata['patient']
    region  = pdata['region']
    tag     = patient+'-'+region

    # process stored data
    
    df_index = pd.read_csv('%s/processed/%s-SU-%s-index.csv' % (HIV_DIR, patient, region), comment='#', memory_map=True)
    df_poly  = pd.read_csv('%s/analysis/%s-poly.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    seqs     = np.loadtxt('%s/%s-poly-seq2state.dat' % (HIV_MPL_DIR, tag), dtype='int')
    times    = np.unique(seqs.T[0])
    
    c_vals   = [C_DEL, C_NEU]
    m_loc    = [[] for s in seqs]
    m_c      = [[] for s in seqs]
    for i in range(2, len(seqs[0])):
        df_temp = df_index[df_index.polymorphic==(i-2)]
        SU  = str(df_temp.iloc[0].SU)
        TF  = str(df_temp.iloc[0].TF)
        loc = int(df_temp.iloc[0].alignment)
        for j in range(len(seqs)):
            n = NUC[seqs[j][i]]
            # SU, and SU does not match TF
            if n==SU and SU!=TF:
                m_loc[j].append(loc)
                m_c[j].append(c_vals[0])
            # Other difference from TF
            elif n!=TF:
                m_loc[j].append(loc)
                m_c[j].append(c_vals[1])

    m_loc = np.array(m_loc)
    m_c   = np.array(m_c)
    
    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN #SLIDE_WIDTH
    hshrink = 1
    goldh   = 0.66 * w * hshrink
    fig     = plt.figure(figsize=(w, goldh))

    box_hist = dict(left=0.15, right=0.85, bottom=0.05, top=0.95)
    gs_hist  = gridspec.GridSpec(len(times), 1, hspace=0.15, wspace=0.15, height_ratios=[np.sum(seqs.T[0]==t) for t in times], **box_hist)
    ax_hist  = [plt.subplot(gs_hist[i, 0]) for i in range(len(times))]

    for k in range(len(times)):
        
        # Subplot - recombination
        
        sp         = 0.2
        rect_props = dict(width=1, height=1, lw=0, alpha=1, ec='none')
        patches    = []
        
        same_time = seqs.T[0]==times[k]
        sub_m_loc = m_loc[same_time]
        sub_m_c   = m_c[same_time]
        
        for i in range(len(sub_m_loc)):
            rect_props['width'] = 1
            for j in range(len(sub_m_loc[i])):
                patches.append(matplotlib.patches.Rectangle(xy=(sub_m_loc[i][j], (1+sp)*i), fc=sub_m_c[i][j], **rect_props))

            start = 0
            for j in range(len(sub_m_loc[i])):
                rect_props['width'] = sub_m_loc[i][j]-start
                patches.append(matplotlib.patches.Rectangle(xy=(start, (1+sp)*i), fc='#f7f7f7', **rect_props))
                start = sub_m_loc[i][j]+1
            rect_props['width'] = 2645-start
            patches.append(matplotlib.patches.Rectangle(xy=(start, (1+sp)*i), fc='#f7f7f7', **rect_props))

        pprops = { 'colors': [BKCOLOR],
                   'xlim': [-1, 2645],
                   'ylim': [ 0, (1+sp)*len(sub_m_loc)],
                   'xticks': [],
                   'yticks': [],
                   'hide': ['top', 'bottom', 'left', 'right']}

        mp.plot(type='scatter', ax=ax_hist[k], x=[[-3]], y=[[-3]], **pprops)

        for patch in patches:
            ax_hist[k].add_artist(patch)

        rec = matplotlib.patches.Rectangle(xy=(0, 0), height=(1+sp)*len(sub_m_loc), width=2645, ec=BKCOLOR, fc='none', lw=AXWIDTH/2, clip_on=False)
        ax_hist[k].add_artist(rec)

    # outside text

    invt = ax_hist[0].transData.inverted()
    xy1  = invt.transform((0,0))
    xy2  = invt.transform((9,9))
    legend_dx = 0.5 * (xy1[0] - xy2[0])
    legend_dy = xy1[1] - xy2[1]
    legend_x  = 2645 + 120
    legend_y  = (1 + 0.2) * np.sum(seqs.T[0]==times[0]) + 0.75 * legend_dy

    rec_props = dict(height=0.5*np.fabs(legend_dy), width=0.5*np.fabs(legend_dx), clip_on=False)
    rec = matplotlib.patches.Rectangle(xy=(legend_x, legend_y + (-0.25*legend_dy)), fc=C_DEL, **rec_props)
    ax_hist[0].add_artist(rec)
    rec = matplotlib.patches.Rectangle(xy=(legend_x, legend_y + ( 1.25*legend_dy)), fc=C_NEU, **rec_props)
    ax_hist[0].add_artist(rec)

    tprops = dict(ha='left', va='center', family=FONTFAMILY, size=SIZELABEL, color=BKCOLOR, clip_on=False)

    ax_hist[0].text(legend_x - legend_dx,  legend_y + (-0.5*legend_dy), 'Superinfecting\nvariants', **tprops)
    ax_hist[0].text(legend_x - legend_dx,  legend_y + ( 1.0*legend_dy), 'Other variants', **tprops)

    tprops['ha'] = 'center'
    for i in range(len(times)):
        ax_hist[i].text(-30, 0.5*(1+sp)*np.sum(seqs.T[0]==times[i]), times[i], rotation=90, **tprops)

    ax_hist[-1].text(0.13,  0.50, 'Time (days)', rotation=90, transform=fig.transFigure, **tprops)
    ax_hist[-1].text(0.50,  0.025, 'HXB2 index', transform=fig.transFigure, **tprops)
    ax_hist[-1].text(   0, -4, '6225', **tprops)
    ax_hist[-1].text(2645, -4, '8794', **tprops)

    plt.savefig('%s/sup-fig-2-cap256-recombination%s' % (FIG_DIR, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('CAP256 recombination done.')


def plot_hive(ax_hive, tag, **pdata):
    """
    Hive plot of \Delta s effects.
    """

    # import stored data
    
    df_ds    = pd.read_csv('%s/analysis/%s-delta-s.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    df_poly  = pd.read_csv('%s/analysis/%s-poly.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    df_index = pd.read_csv('%s/processed/%s-index.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    L        = len(df_index)
    
    # set epitope labels

    epitope_start = []
    epitope_end   = []
    epitope_label = list(set(list(df_index[pd.notnull(df_index.epitope)].epitope)))
    print(tag, epitope_label)
    
    for e in epitope_label:
        epitope_start.append(np.min(df_index[df_index.epitope==e].alignment))
        epitope_end.append(np.max(df_index[df_index.epitope==e].alignment))

    # PLOT FIGURE

    # populate links

    idx_mask_pos  = []
    idx_delta_pos = []
    ds_pos        = []
    idx_mask_neg  = []
    idx_delta_neg = []
    ds_neg        = []

    c_cutoff  = 0.003
    ds_cutoff = 0.004
    ds_max    = 0.010

    df_ds = df_ds[np.fabs(df_ds.effect)>ds_cutoff]
    for it, entry in df_ds.iterrows():
        if (entry.mask_polymorphic_index==entry.target_polymorphic_index#):
#            continue
            and entry.mask_nucleotide==entry.target_nucleotide):
            continue
        if entry.effect>0:
            mask_alignment_index   = df_index[df_index.polymorphic==  entry.mask_polymorphic_index].iloc[0].alignment
            target_alignment_index = df_index[df_index.polymorphic==entry.target_polymorphic_index].iloc[0].alignment
            idx_mask_pos.append(mask_alignment_index)
            idx_delta_pos.append(target_alignment_index)
            ds_pos.append(entry.effect)
        elif entry.effect<0:
            mask_alignment_index   = df_index[df_index.polymorphic==  entry.mask_polymorphic_index].iloc[0].alignment
            target_alignment_index = df_index[df_index.polymorphic==entry.target_polymorphic_index].iloc[0].alignment
            idx_mask_neg.append(mask_alignment_index)
            idx_delta_neg.append(target_alignment_index)
            ds_neg.append(entry.effect)
    
    ds_pos           = (np.array(ds_pos) - c_cutoff) / ds_max
    ds_pos[ds_pos>1] = 1
    pos_sort         = np.argsort(ds_pos)[::-1]
    idx_mask_pos     = np.array(idx_mask_pos)[pos_sort]
    idx_delta_pos    = np.array(idx_delta_pos)[pos_sort]
    ds_pos           = ds_pos[pos_sort]

    ds_neg           = (np.fabs(np.array(ds_neg)) - c_cutoff) / ds_max
    ds_neg[ds_neg>1] = 1
    neg_sort         = np.argsort(ds_neg)[::-1]
    idx_mask_neg     = np.array(idx_mask_neg)[neg_sort]
    idx_delta_neg    = np.array(idx_delta_neg)[neg_sort]
    ds_neg           = ds_neg[neg_sort]

    # arc plot for large values of Delta s

    c_pos  = LCOLOR
    c_neg  = LCOLOR
    c_circ = { True : c_neg, False : c_pos }

    r_min     = 0.05
    r_max     = 1.00
    r_norm    = float(L) / (r_max - r_min)
    arc_mult  = SIZELINE * ds_max / 1.5
    arc_alpha = 1

    mask_angle_p =  -np.pi/2
    mask_angle_n = 3*np.pi/2
    ds_pos_angle =   np.pi/6
    ds_neg_angle = 5*np.pi/6
    da           = 0.005

    circ_color  = [hls_to_rgb(0.02, 0.53 * ds + 1. * (1 - ds), 0.83) for ds in ds_pos]
    circ_color += [hls_to_rgb(0.58, 0.53 * ds + 1. * (1 - ds), 0.60) for ds in ds_neg]
    circ_rad    = [[r_min + idx_mask_pos[i]/r_norm, r_min + idx_delta_pos[i]/r_norm] for i in range(len(ds_pos))]
    circ_rad   += [[r_min + idx_mask_neg[i]/r_norm, r_min + idx_delta_neg[i]/r_norm] for i in range(len(ds_neg))]
    circ_arc    = [dict(lw=arc_mult * np.fabs(0.1) / ds_cutoff, alpha=arc_alpha) for ds in ds_pos]
    circ_arc   += [dict(lw=arc_mult * np.fabs(0.1) / ds_cutoff, alpha=arc_alpha) for ds in ds_neg]
    circ_angle  = [[mask_angle_p+da/r[0], ds_pos_angle-da/r[1]] for r in circ_rad[:len(ds_pos)]]
    circ_angle += [[mask_angle_n-da/r[0], ds_neg_angle+da/r[1]] for r in circ_rad[len(ds_pos):]]
    circ_bezrad = [(r[0]+r[1])/2 for r in circ_rad]
    circ_x      = [i for i in ds_pos] + [i for i in ds_neg]
    circ_y      = [i for i in ds_pos] + [i for i in ds_neg]
    
    # mark epitopes

    for i in range(len(epitope_label)):
        r_mid = r_min + (epitope_start[i] + epitope_end[i]) / (2 * r_norm)
        scatter_x = [r_mid * np.cos(mask_angle_p),
                     r_mid * np.cos(ds_pos_angle),
                     r_mid * np.cos(ds_neg_angle) ]
        scatter_y = [r_mid * np.sin(mask_angle_p),
                     r_mid * np.sin(ds_pos_angle),
                     r_mid * np.sin(ds_neg_angle) ]
        smallprops = dict(lw=0, marker='o', s=1.0*SMALLSIZEDOT, zorder=9999, clip_on=False)
        bigprops   = dict(lw=0, marker='o', s=1.7*SMALLSIZEDOT, zorder=9998, clip_on=False)
        mp.scatter(ax=ax_hive, x=[scatter_x], y=[scatter_y], colors=['#FFFFFF'], plotprops=smallprops)
        mp.scatter(ax=ax_hive, x=[scatter_x], y=[scatter_y], colors=[BKCOLOR],   plotprops=  bigprops)

    # Make plot

    pprops = { 'colors':   circ_color,
               'xlim':     [-1.0, 1.0],
               'ylim':     [-1.0, 1.0],
               'size':     L,
               'rad':      circ_rad,
               'arcprops': circ_arc,
               'angle':    circ_angle,
               'bezrad':   circ_bezrad,
               'noaxes':   True }

    plotprops = dict(lw=AXWIDTH, ls='-', zorder=999)
    line_x = [[r_min * np.cos(mask_angle_p), r_max * np.cos(mask_angle_p)],
              [r_min * np.cos(ds_pos_angle), r_max * np.cos(ds_pos_angle)],
              [r_min * np.cos(ds_neg_angle), r_max * np.cos(ds_neg_angle)] ]
    line_y = [[r_min * np.sin(mask_angle_p), r_max * np.sin(mask_angle_p)],
              [r_min * np.sin(ds_pos_angle), r_max * np.sin(ds_pos_angle)],
              [r_min * np.sin(ds_neg_angle), r_max * np.sin(ds_neg_angle)] ]
    mp.line(ax=ax_hive, x=line_x, y=line_y, colors=[BKCOLOR for i in range(len(line_x))], plotprops=plotprops)

    mp.plot(type='circos', ax=ax_hive, x=circ_x, y=circ_y, **pprops)
    
    
def plot_supplementary_figure_s_conditions(**pdata):
    """
    Scatter plot of selection coefficients inferred using different conditions for
    data processing.
    """
    
    # Read in data

    columns = pdata['columns']
    labels  = pdata['labels']
    
    df = pd.read_csv('%s/analysis/total-selection-merged.csv' % (HIV_DIR), comment='#', memory_map=True)

    # Define plot

    w     = DOUBLE_COLUMN
    goldh = w

    fig = plt.figure(figsize=(w, goldh))
    box = dict(left=0.08, right=0.99, bottom=0.08, top=0.99)
    gs  = gridspec.GridSpec(len(columns)-1, len(columns)-1, **box)
    ax  = [[0 for j in range(len(columns)-1)] for i in range(len(columns)-1)]
        
    for i in range(len(columns)):
        for j in range(i+1,len(columns)):
            ax[i][j-1] = plt.subplot(gs[j-1, i])

    # MAKE SUBPLOTS

    scatterprops = dict(lw=0, s=SMALLSIZEDOT*1.0, marker='o', alpha=0.2, c=C_MPL)

    pprops = { 'xlim':  [-0.10, 0.15],
               'ylim':  [-0.10, 0.15],
               'xticks': [],
               'yticks': [],
               'theme': 'open',
               #'hide':  ['left', 'bottom'],
               'plotprops': scatterprops }

    print('choice 1\t\tchoice 2\t\tR\tR^2\tp')
    
    pearson_r2s = []
    spearman_rs = []
    diffs = []
    abs_s = []

    for i in range(len(columns)):
        for j in range(i+1,len(columns)):
        
            tempprops = deepcopy(pprops)
            if j==len(columns)-1:
                tempprops['xlabel'] = labels[i]
                tempprops['xticks'] = [-0.10, 0, 0.10]
                tempprops['xticklabels'] = [-10, 0, 10]
            if i==0:
                tempprops['ylabel'] = labels[j]
                tempprops['yticks'] = [-0.10, 0, 0.10]
                tempprops['yticklabels'] = [-10, 0, 10]
            
            df_sub = df[['variant', columns[j], columns[i]]].dropna()
            x = df_sub[columns[j]]
            y = df_sub[columns[i]]
            diffs += list(np.fabs(x - y))
            abs_s += list(np.fabs(x)) + list(np.fabs(y))
            
#            # Identify large differences
#            df_diff = df_sub[np.fabs(df_sub[columns[j]]-df_sub[columns[i]])>0.03]
#            for df_iter, df_entry in df_diff.iterrows():
#                print('\t%s\t%s %.2f\t%s %.2f' % (df_entry.variant, columns[i], df_entry[columns[i]], columns[j], df_entry[columns[j]]))
            
            r_value, p_value = st.pearsonr(x, y)
            pearson_r2s.append(r_value**2)
            print('%s\t\t%s\t\t%.2f\t%.2f\t%.2e' % (columns[i], columns[j], r_value, r_value**2, p_value))
            
            r_value, p_value = st.spearmanr(x, y)
            spearman_rs.append(r_value)
            
            mp.plot(type='scatter', ax=ax[i][j-1], x=[x], y=[y], **tempprops)
            
    # Aggregate correlation statistics
    
    print('')
    print('metric\t\tminimum\t\tmaximum\t\tmean\t\tmedian')
    print('Pearson R2\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f' % (np.min(pearson_r2s), np.max(pearson_r2s), np.mean(pearson_r2s), np.median(pearson_r2s)))
    print('Spearman r\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f' % (np.min(spearman_rs), np.max(spearman_rs), np.mean(spearman_rs), np.median(spearman_rs)))
    print('|Delta s| \t%.2f\t\t%.2f\t\t%.2e\t%.2e' % (np.min(diffs), np.max(diffs), np.mean(diffs), np.median(diffs)))
    print('|s|       \t%.2f\t\t%.2f\t\t%.2e\t%.2e' % (np.min(abs_s), np.max(abs_s), np.mean(abs_s), np.median(abs_s)))
            
    # Save data

    plt.savefig('%s/ed-fig-10-s-conditions%s' % (FIG_DIR, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('\nSelection coefficient comparison done.')


def plot_supplementary_figure_varying_N(**pdata):
    """
    Example evolutionary trajectory for a 50-site system and inferred selection coefficients,
    together with aggregate properties across sampling levels.
    """
    
    # unpack passed data

    n_gen  = pdata['n_gen']
    dg     = pdata['dg']
    N      = pdata['N']
    xfile  = pdata['xfile']
    method = pdata['method']

    n_ben = pdata['n_ben']
    n_neu = pdata['n_neu']
    n_del = pdata['n_del']
    s_ben = pdata['s_ben']
    s_neu = pdata['s_neu']
    s_del = pdata['s_del']

    # load and process data files

    data  = np.loadtxt('%s/data/%s.dat' % (WFS_DIR, xfile))
    times = np.unique(data.T[0])
    x     = []
    for i in range(0, n_gen, dg):
        idx    = data.T[0]==i
        t_data = data[idx].T[2:].T
        t_num  = data[idx].T[1].T
        t_freq = np.einsum('i,ij->j', t_num, t_data) / float(np.sum(t_num))
        x.append(t_freq)
    x = np.array(x).T

    s_true = [s_ben for i in range(n_ben)] + [0 for i in range(n_neu)] + [s_del for i in range(n_del)]
    s_inf  = np.loadtxt('%s/%s_%s.dat' % (SIM_MPL_DIR, xfile.split('wfsim_')[1], method))
    cov    = np.loadtxt('%s/covariance-%s.dat' % (SIM_MPL_DIR, xfile.split('wfsim_')[1]))
    ds     = np.linalg.inv(cov) / N

    # PLOT FIGURE

    ## set up figure grid

    w     = ONE_FIVE_COLUMN
    goldh = w * 1.5
    fig   = plt.figure(figsize=(w, goldh))

    n_rows = [   n_ben,    n_neu,       n_del]
    offset = [       0,    n_ben, n_ben+n_neu]
    tag    = [   'ben',    'neu',       'del']
    colors = [   C_BEN,    C_NEU,       C_DEL]
    fc     = [C_BEN_LT, C_NEU_LT,    C_DEL_LT]

    htot = 0.65 + 0.01
    dh   = (htot - 0.08) / float(n_ben + n_neu + n_del + 2)

    boxl = [0.12, 0.12, 0.12, 0.12]
    boxr = [0.61, 0.61, 0.61, 0.61]
    boxb = [htot - (dh * n_ben), htot - (dh * (n_ben + n_neu + 1)), htot - (dh * (n_ben + n_neu + n_del + 2))]
    boxt = [               htot,         htot - (dh * (n_ben + 1)), htot - (dh * (n_ben + n_neu + 2))]
    gs   = [0 for k in range(len(tag))]
    
    ## a -- population size vs. time
    
    top = 0.97
    dx  = 0.03
    dy  = 0.05
    tempgs = gridspec.GridSpec(1, 1)
    tempgs.update(left=boxl[0], right=top, bottom=top-2*dy, top=top)
    ax     = plt.subplot(tempgs[0, 0])

    lineprops = { 'lw' : SIZELINE*1.5, 'ls' : '-', 'alpha' : 1 }

    pprops = { 'xticks'      : [0, 100, 200, 300, 400, 500],
               'yticks'      : [1e0, 1e3, 1e6],
               'logy'        : True,
               'xlabel'      : 'Generation',
               'ylabel'      : 'Population\nsize, '+r'$N$',
               'plotprops'   : lineprops,
               'axoffset'    : 0.1,
               'theme'       : 'open' }
               
    def get_N(t):
        """ Return population size as a function of time. """
        
        N0 = 1
        N1 = 1000000
        N2 = 1000
        t0 = 0
        t1 = 50
        t2 = 100
        e_rise = np.log(N1)/(t1-t0)
        e_fall = np.log(N2/N1)/(t2-t1)
        if t==t0:
            return N0
        elif t>t0 and t<=t1:
            return np.round(N0 * np.exp(e_rise * (t - t0)))
        elif t>t1 and t<=t2:
            return np.round(N1 * np.exp(e_fall * (t - t1)))
        else:
            return N2

    xdat = [np.arange(0, 500, 1)]
    ydat = [[get_N(t) for t in xdat[0]]]
    mp.plot(type='line', ax=ax, x=xdat, y=ydat, colors=[C_MPL], **pprops)
    
    ax.text(boxl[0]-dx-0.03, top+0.01, 'a'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## b -- all trajectories together

    dx = 0.03
    dy = 0.05
    tempgs = gridspec.GridSpec(1, 1)
    tempgs.update(left=boxl[0], right=boxr[0], bottom=top-5*dy, top=top-3*dy)
    ax     = plt.subplot(tempgs[0, 0])

    lineprops = { 'lw' : SIZELINE*1.5, 'ls' : '-', 'alpha' : 0.6 }

    pprops = { 'xticks'      : [0, 100, 200, 300, 400, 500],
               'yticks'      : [0, 1],
               'yminorticks' : [0.25, 0.5, 0.75],
               'nudgey'      : 1.1,
               'xlabel'      : 'Generation',
               'ylabel'      : 'Allele\nfrequency, '+r'$x$',
               'plotprops'   : lineprops,
               'axoffset'    : 0.1,
               'theme'       : 'open' }

    xdat = [range(0, n_gen, dg) for k in range(len(x))]
    ydat = [k for k in x]
    mp.plot(type='line', ax=ax, x=xdat, y=ydat, colors=[C_BEN_LT for k in range(n_ben)] + [LCOLOR for k in range(n_neu)] + [C_DEL_LT for k in range(n_del)], **pprops)
#    mp.plot(type='line', ax=ax, x=xdat, y=ydat, colors=[LCOLOR for k in range(len(x))], **pprops)

    ## c -- individual beneficial/neutral/deleterious trajectories and selection coefficients

    idx = 0
    for k in range(len(tag)):
        
        ### trajectories
        
        gs[k] = gridspec.GridSpec(n_rows[k], 2)
        gs[k].update(left=boxl[k], right=boxr[k], bottom=boxb[k], top=boxt[k], wspace=0.05)
        ax = [[plt.subplot(gs[k][i, 0]), plt.subplot(gs[k][i, 1])] for i in range(n_rows[k])]
        
        legendprops = { 'loc' : 4, 'frameon' : False, 'scatterpoints' : 1, 'handletextpad' : 0.1,
                        'prop' : {'size' : SIZELABEL}, 'ncol' : 1 }
        lineprops   = { 'lw' : SIZELINE*1.5, 'linestyle' : '-', 'alpha' : 1.0 }
        fillprops   = { 'lw' : 0, 'alpha' : 0.3, 'interpolate' : True }
        
        pprops = { 'xticks' : [],
                   'yticks' : [],
                   'hide'   : ['top','bottom','left','right'],
                   'theme'  : 'open' }
        
        for i in range(n_rows[k]):
            pprops['colors'] = [colors[k]]
            pprops['xlim']   = [    0,  500]
            pprops['ylim']   = [-0.08, 1.08]
            ydat             = x[offset[k]+i]
            if (i==n_rows[k]-1) and (k==len(tag)-1):
                pprops['xticks']   = [0, 100, 200, 300, 400, 500]
                pprops['xlabel']   = 'Generation'
                pprops['hide']     = ['top','left','right']
                pprops['axoffset'] = 0.3
            mp.line(             ax=ax[i][0], x=[range(0, n_gen, dg)], y=[ydat], plotprops=lineprops, **pprops)
            mp.plot(type='fill', ax=ax[i][0], x=[range(0, n_gen, dg)], y=[ydat], plotprops=fillprops, **pprops)
        
        ### selection coefficient estimates
        
        sprops = { 'lw' : 0, 's' : 9., 'marker' : 'o' }
        
        pprops = { 'yticks'  : [],
            'xticks'  : [],
            'hide'    : ['top','bottom','left','right'],
            'theme'   : 'open' }
        
        for i in range(n_rows[k]):
            pprops['xlim'] = [-0.04, 0.04]
            pprops['ylim'] = [  0.5,  1.5]
            ydat           = [1]
            xdat           = [s_inf[offset[k]+i]]
            xerr           = np.sqrt(ds[offset[k]+i][offset[k]+i])
            if (i==n_rows[k]-1) and (k==len(tag)-1):
                pprops['xticks']      = [-0.04,   -0.02,      0,   0.02,   0.04]
                pprops['xticklabels'] = [  ' ', r'$-2$', r'$0$', r'$2$', r'$4$']
                pprops['xlabel']      = 'Inferred selection\ncoefficient, ' + r'$\hat{s}$' + ' (%)'
                pprops['hide']        = ['top','left','right']
                pprops['axoffset']    = 0.3
            mp.plot(type='error', ax=ax[i][1], x=[xdat], y=[ydat], xerr=[xerr], colors=[colors[k]], **pprops)
            ax[i][1].axvline(x=s_true[idx], ls=':', lw=SIZELINE, color=BKCOLOR)
            idx += 1

    ### bounding boxes

    ax[0][0].text(boxl[0]-0.06, top-3*dy+0.01, 'b'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax[0][0].text(boxl[0]-0.06, boxt[0], 'c'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    lineprops = { 'lw' : AXWIDTH/2., 'ls' : '-', 'alpha' : 1.0 }
    pprops    = { 'xlim' : [0, 1], 'xticks' : [], 'ylim' : [0, 1], 'yticks' : [],
        'hide' : ['top','bottom','left','right'], 'plotprops' : lineprops }
    txtprops = {'ha' : 'right', 'va' : 'center', 'color' : BKCOLOR, 'family' : FONTFAMILY,
        'size' : SIZELABEL, 'rotation' : 90, 'transform' : fig.transFigure}
    ax[0][0].text(boxl[0]-0.01, (boxb[0]+boxt[0])/2.,  'Beneficial', **txtprops)
    ax[0][0].text(boxl[0]-0.01, (boxb[1]+boxt[1])/2.,     'Neutral', **txtprops)
    ax[0][0].text(boxl[0]-0.01, (boxb[2]+boxt[2])/2., 'Deleterious', **txtprops)

    boxprops = {'ec' : BKCOLOR, 'lw' : SIZELINE/2., 'fc' : 'none', 'clip_on' : False, 'zorder' : -100}

    dx  = 0.005
    dxl = 0.000
    dxr = 0.001
    dy  = 0.003
    for k in range(len(tag)):
        ll = boxl[k] + dxl                     # left box left
        lb = boxb[k] - dy                      # left box bottom
        rl = (boxl[k] + boxr[k])/2. + dx + dxl # right box left
        wd = (boxr[k] - boxl[k])/2. - dxr      # box width
        ht = (boxt[k] - boxb[k]) + (2. * dy)   # box height
        
        recL = matplotlib.patches.Rectangle(xy=(ll, lb), width=wd, height=ht, transform=fig.transFigure, **boxprops)
        recL = ax[0][0].add_patch(recL)
        recR = matplotlib.patches.Rectangle(xy=(rl, lb), width=wd, height=ht, transform=fig.transFigure, **boxprops)
        recR = ax[0][0].add_patch(recR)

    ## d -- histogram of selection coefficients

    ### set up grid

    box_l  = dict(left=0.72, right=top, bottom=0.05, top=top-3*0.05-0.01)
    #box_ru = dict(left=0.69, right=0.97, bottom=0.63, top=0.96)
    #box_rl = dict(left=0.69, right=0.97, bottom=0.14, top=0.47)

    gs_l  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_l)
    #gs_ru = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_ru)
    #gs_rl = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_rl)
    ax_l  = plt.subplot(gs_l[0, 0])
    #ax_ru = plt.subplot(gs_ru[0, 0])
    #ax_rl = plt.subplot(gs_rl[0, 0])

    ### plot histogram

    df   = pd.read_csv('%s/MPL_simple_N_collected_extended.csv.gz' % (SIM_DIR), memory_map=True)
    df   = df[df.method==method]
    df_s = df[(df.deltat==10) & (df.ns==100)]

    ben_cols = ['s%d' % i for i in range(n_ben)]
    neu_cols = ['s%d' % i for i in range(n_ben, n_ben+n_neu)]
    del_cols = ['s%d' % i for i in range(n_ben+n_neu, n_ben+n_neu+n_del)]

    colors     = [C_BEN, C_NEU, C_DEL]
    tags       = ['beneficial', 'neutral', 'deleterious']
    cols       = [ben_cols, neu_cols, del_cols]
    s_true_loc = [s_ben, s_neu, s_del]

    dashlineprops = { 'lw' : SIZELINE * 2.0, 'ls' : ':', 'alpha' : 0.5, 'color' : BKCOLOR }
    histprops = dict(histtype='bar', lw=SIZELINE/2, rwidth=0.8, ls='solid', alpha=0.7, edgecolor='none',
                     orientation='horizontal')
    pprops = { 'xlim'        : [-0.04, 0.04],
               'xticks'      : [  -0.04,   -0.03,   -0.02,   -0.01,     0.,   0.01,   0.02,   0.03,   0.04],
               'yticklabels' : [r'$-4$', r'$-3$', r'$-2$', r'$-1$', r'$0$', r'$1$', r'$2$', r'$3$', r'$4$'],
               'ylim'        : [0., 0.075],
               'yticks'      : [0., 0.025, 0.05, 0.075],
               'xlabel'      : 'Frequency',
               'ylabel'      : 'Inferred selection coefficient, ' + r'$\hat{s}$' + ' (%)',
               'bins'        : np.arange(-0.04, 0.04+0.001, 0.001),
               'combine'     : True,
               'plotprops'   : histprops,
               'axoffset'    : 0.1,
               'theme'       : 'boxed' }

    for i in range(len(tags)):
        x = [np.array(df_s[cols[i]]).flatten()]
        tprops = dict(ha='left', va='center', family=FONTFAMILY, size=SIZELABEL, rotation=270, clip_on=False)
        ax_l.text(0.102, s_true_loc[i], r'$s_{%s}$' % (tags[i]), color=colors[i], **tprops)
        ax_l.axhline(y=s_true_loc[i], **dashlineprops)
        if i<len(tags)-1: mp.hist(             ax=ax_l, x=x, colors=[colors[i]], **pprops)
        else:             mp.plot(type='hist', ax=ax_l, x=x, colors=[colors[i]], **pprops)

    ## outside text labels

    tprops = dict(color=BKCOLOR, ha='center', va='center', family=FONTFAMILY, size=SIZELABEL,
                  clip_on=False, transform=fig.transFigure)
    dx = -0.03
    dy =  0.02

    ax_l.text(box_l['left']+dx,  box_l['top']+dy, 'd'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/revision/varying-N%s' % (FIG_DIR, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('MPL supplementary example (varying N) done.')


def plot_supplementary_figure_varying_s(**pdata):
    """
    Example evolutionary trajectory for a 50-site system and inferred selection coefficients,
    together with aggregate properties across sampling levels.
    """
    
    # unpack passed data

    n_gen  = pdata['n_gen']
    dg     = pdata['dg']
    N      = pdata['N']
    method = pdata['method']
    name   = pdata['name']

    n_ben = pdata['n_ben']
    n_neu = pdata['n_neu']
    n_del = pdata['n_del']
    s_ben = pdata['s_ben']
    s_neu = pdata['s_neu']
    s_del = pdata['s_del']

    # load and process data files

    s_true = [s_ben for i in range(n_ben-1)] + [0 for i in range(n_neu)] + [s_del for i in range(n_del)]
    s_inf  = np.load('%s/single-%s.npz' % (SIM_DIR, name))['selection_constant']
    cov    = np.load('%s/WF-single-%s.npz' % (SIM_DIR, name))['covar_int'] / N
    ds     = np.linalg.inv(cov) / N

    # PLOT FIGURE

    ## set up figure grid

    w     = ONE_FIVE_COLUMN
    goldh = w * 1.5
    fig   = plt.figure(figsize=(w, goldh))

    n_rows = [ n_ben-1,       n_neu,    n_del]
    offset = [       1, n_ben+n_del,    n_ben]
    tag    = [   'ben',       'neu',    'del']
    colors = [   C_BEN,       C_NEU,    C_DEL]
    fc     = [C_BEN_LT,    C_NEU_LT, C_DEL_LT]

    htot = 0.65 + 0.01
    dh   = (htot - 0.08) / float(n_ben + n_neu + n_del + 2)

    boxl = [0.12, 0.12, 0.12, 0.12]
    boxr = [0.61, 0.61, 0.61, 0.61]
    boxb = [htot - (dh * n_ben), htot - (dh * (n_ben + n_neu + 1)), htot - (dh * (n_ben + n_neu + n_del + 2))]
    boxt = [               htot,         htot - (dh * (n_ben + 1)), htot - (dh * (n_ben + n_neu + 2))]
    gs   = [0 for k in range(len(tag))]
    
    ## a -- selection coefficient vs. time
    
    top = 0.97
    dx  = 0.03
    dy  = 0.05
    tempgs = gridspec.GridSpec(1, 1)
    tempgs.update(left=boxl[0], right=boxr[0], bottom=top-2*dy, top=top)
    ax     = plt.subplot(tempgs[0, 0])

    lineprops = { 'lw' : SIZELINE*1.5, 'ls' : '-', 'alpha' : 1 }

    pprops = { 'xticks'      : [0, 200, 400, 600, 800, 1000],
               'yticks'      : [0, 0.02, 0.04],
               'yticklabels' : [r'$0$', r'$2$', r'$4$'],
               'xlabel'      : 'Generation',
               'ylabel'      : 'Time-varying selection\ncoefficient, '+r'$s$' + ' (%)',
               'plotprops'   : lineprops,
               'axoffset'    : 0.1,
               'theme'       : 'open' }
               
    xdat = [np.arange(0, 1000, 1)]
    ydat = [np.load('%s/s-single-%s.npy' % (SIM_DIR, name)).T[0]]
    savg = np.mean(ydat[0]) if ydat[0][0]==0 else ydat[0][0]
    mp.plot(type='line', ax=ax, x=xdat, y=ydat, colors=[C_MPL], **pprops)
    ax.axhline(y=savg, ls=':', lw=SIZELINE, color=BKCOLOR)
    
    ax.text(boxl[0]-dx-0.03, top+0.01, 'a'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    
    ## b -- time-varying selection coefficients
    
    box_t  = dict(left=0.72, right=top, bottom=top-2*0.05, top=top-0.01)
    gs_t  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_t)
    ax_t  = plt.subplot(gs_t[0, 0])

    ### plot histogram

    df   = pd.read_csv('%s/MPL_%s_collected_extended.csv.gz' % (SIM_DIR, name), memory_map=True)
    df   = df[df.method==method]
    df_s = df[(df.deltat==10) & (df.ns==100)]

    histprops = dict(histtype='bar', lw=SIZELINE/2, rwidth=0.8, ls='solid', alpha=0.7, edgecolor='none',
                     orientation='horizontal')
    pprops = { 'xlim'        : [0, 0.04],
               'xticks'      : [0, 0.02, 0.04],
               'yticklabels' : [r'$0$', r'$2$', r'$4$'],
               'ylim'        : [0., 0.075],
               'yticks'      : [0., 0.025, 0.05, 0.075],
               'xlabel'      : 'Frequency',
               'ylabel'      : 'Inferred selection\ncoefficient, ' + r'$\hat{s}$' + ' (%)',
               'bins'        : np.arange(0, 0.04+0.001, 0.001),
               'combine'     : True,
               'plotprops'   : histprops,
               'axoffset'    : 0.1,
               'theme'       : 'boxed' }

    x = [np.array(df_s['s0']).flatten()]
    mp.plot(type='hist', ax=ax_t, x=x, colors=[C_MPL], **pprops)
    ax_t.axhline(y=savg, ls=':', lw=SIZELINE, color=BKCOLOR)
    
    ax.text(0.68, top+0.01, 'b'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## c -- all trajectories together

    dx = 0.03
    dy = 0.05
    tempgs = gridspec.GridSpec(1, 1)
    tempgs.update(left=boxl[0], right=boxr[0], bottom=top-5*dy, top=top-3*dy)
    ax     = plt.subplot(tempgs[0, 0])

    lineprops = { 'lw' : SIZELINE*1.5, 'ls' : '-', 'alpha' : 0.6 }

    pprops = { 'xticks'      : [0, 200, 400, 600, 800, 1000],
               'yticks'      : [0, 1],
               'yminorticks' : [0.25, 0.5, 0.75],
               'nudgey'      : 1.1,
               'xlabel'      : 'Generation',
               'ylabel'      : 'Allele\nfrequency, '+r'$x$',
               'plotprops'   : lineprops,
               'axoffset'    : 0.1,
               'theme'       : 'open' }
               
    data = np.load('%s/WF-single-%s.npz' % (SIM_DIR, name))['traj']
    x    = data.T

    xdat = [range(0, n_gen, dg) for k in range(len(x))]
    ydat = [k for k in x]
    mp.plot(type='line', ax=ax, x=xdat, y=ydat,
            colors=[C_MPL] + [C_BEN_LT for k in range(n_ben-1)] + [C_DEL_LT for k in range(n_del)] + [LCOLOR for k in range(n_neu)], **pprops)
#    mp.plot(type='line', ax=ax, x=xdat, y=ydat, colors=[LCOLOR for k in range(len(x))], **pprops)

    ## d -- individual beneficial/neutral/deleterious trajectories and selection coefficients

    idx = 0
    for k in range(len(tag)):
        
        ### trajectories
        
        gs[k] = gridspec.GridSpec(n_rows[k], 2)
        gs[k].update(left=boxl[k], right=boxr[k], bottom=boxb[k], top=boxt[k], wspace=0.05)
        ax = [[plt.subplot(gs[k][i, 0]), plt.subplot(gs[k][i, 1])] for i in range(n_rows[k])]
        
        legendprops = { 'loc' : 4, 'frameon' : False, 'scatterpoints' : 1, 'handletextpad' : 0.1,
                        'prop' : {'size' : SIZELABEL}, 'ncol' : 1 }
        lineprops   = { 'lw' : SIZELINE*1.5, 'linestyle' : '-', 'alpha' : 1.0 }
        fillprops   = { 'lw' : 0, 'alpha' : 0.3, 'interpolate' : True }
        
        pprops = { 'xticks' : [],
                   'yticks' : [],
                   'hide'   : ['top','bottom','left','right'],
                   'theme'  : 'open' }
        
        for i in range(n_rows[k]):
            pprops['colors'] = [colors[k]]
            pprops['xlim']   = [    0, 1000]
            pprops['ylim']   = [-0.08, 1.08]
            ydat             = x[offset[k]+i]
            if (i==n_rows[k]-1) and (k==len(tag)-1):
                pprops['xticks']   = [0, 200, 400, 600, 800, 1000]
                pprops['xlabel']   = 'Generation'
                pprops['hide']     = ['top','left','right']
                pprops['axoffset'] = 0.3
            mp.line(             ax=ax[i][0], x=[range(0, n_gen, dg)], y=[ydat], plotprops=lineprops, **pprops)
            mp.plot(type='fill', ax=ax[i][0], x=[range(0, n_gen, dg)], y=[ydat], plotprops=fillprops, **pprops)
        
        ### selection coefficient estimates
        
        sprops = { 'lw' : 0, 's' : 9., 'marker' : 'o' }
        
        pprops = { 'yticks'  : [],
            'xticks'  : [],
            'hide'    : ['top','bottom','left','right'],
            'theme'   : 'open' }
        
        for i in range(n_rows[k]):
            pprops['xlim'] = [-0.04, 0.04]
            pprops['ylim'] = [  0.5,  1.5]
            ydat           = [1]
            xdat           = [s_inf[offset[k]+i]]
            xerr           = np.sqrt(ds[offset[k]+i][offset[k]+i])
            if (i==n_rows[k]-1) and (k==len(tag)-1):
                pprops['xticks']      = [-0.04,   -0.02,      0,   0.02,   0.04]
                pprops['xticklabels'] = [  ' ', r'$-2$', r'$0$', r'$2$', r'$4$']
                pprops['xlabel']      = 'Inferred selection\ncoefficient, ' + r'$\hat{s}$' + ' (%)'
                pprops['hide']        = ['top','left','right']
                pprops['axoffset']    = 0.3
            mp.plot(type='error', ax=ax[i][1], x=[xdat], y=[ydat], xerr=[xerr], colors=[colors[k]], **pprops)
            ax[i][1].axvline(x=s_true[idx], ls=':', lw=SIZELINE, color=BKCOLOR)
            idx += 1

    ### bounding boxes

    ax[0][0].text(boxl[0]-0.06, top-3*dy+0.01, 'c'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax[0][0].text(boxl[0]-0.06, boxt[0], 'd'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    lineprops = { 'lw' : AXWIDTH/2., 'ls' : '-', 'alpha' : 1.0 }
    pprops    = { 'xlim' : [0, 1], 'xticks' : [], 'ylim' : [0, 1], 'yticks' : [],
        'hide' : ['top','bottom','left','right'], 'plotprops' : lineprops }
    txtprops = {'ha' : 'right', 'va' : 'center', 'color' : BKCOLOR, 'family' : FONTFAMILY,
        'size' : SIZELABEL, 'rotation' : 90, 'transform' : fig.transFigure}
    ax[0][0].text(boxl[0]-0.01, (boxb[0]+boxt[0])/2.,  'Beneficial', **txtprops)
    ax[0][0].text(boxl[0]-0.01, (boxb[1]+boxt[1])/2.,     'Neutral', **txtprops)
    ax[0][0].text(boxl[0]-0.01, (boxb[2]+boxt[2])/2., 'Deleterious', **txtprops)

    boxprops = {'ec' : BKCOLOR, 'lw' : SIZELINE/2., 'fc' : 'none', 'clip_on' : False, 'zorder' : -100}

    dx  = 0.005
    dxl = 0.000
    dxr = 0.001
    dy  = 0.003
    for k in range(len(tag)):
        ll = boxl[k] + dxl                     # left box left
        lb = boxb[k] - dy                      # left box bottom
        rl = (boxl[k] + boxr[k])/2. + dx + dxl # right box left
        wd = (boxr[k] - boxl[k])/2. - dxr      # box width
        ht = (boxt[k] - boxb[k]) + (2. * dy)   # box height
        
        recL = matplotlib.patches.Rectangle(xy=(ll, lb), width=wd, height=ht, transform=fig.transFigure, **boxprops)
        recL = ax[0][0].add_patch(recL)
        recR = matplotlib.patches.Rectangle(xy=(rl, lb), width=wd, height=ht, transform=fig.transFigure, **boxprops)
        recR = ax[0][0].add_patch(recR)

    ## e -- histogram of selection coefficients

    ### set up grid

    box_l  = dict(left=0.72, right=top, bottom=0.05, top=top-3*0.05-0.01)
    #box_ru = dict(left=0.69, right=0.97, bottom=0.63, top=0.96)
    #box_rl = dict(left=0.69, right=0.97, bottom=0.14, top=0.47)

    gs_l  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_l)
    #gs_ru = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_ru)
    #gs_rl = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_rl)
    ax_l  = plt.subplot(gs_l[0, 0])
    #ax_ru = plt.subplot(gs_ru[0, 0])
    #ax_rl = plt.subplot(gs_rl[0, 0])

    ### plot histogram

    df   = pd.read_csv('%s/MPL_%s_collected_extended.csv.gz' % (SIM_DIR, name), memory_map=True)
    df   = df[df.method==method]
    df_s = df[(df.deltat==10) & (df.ns==100)]

    ben_cols = ['s%d' % i for i in range(1, n_ben)]
    del_cols = ['s%d' % i for i in range(n_ben, n_ben+n_del)]
    neu_cols = ['s%d' % i for i in range(n_ben+n_del, n_ben+n_del+n_neu)]

    colors     = [C_BEN, C_NEU, C_DEL]
    tags       = ['beneficial', 'neutral', 'deleterious']
    cols       = [ben_cols, neu_cols, del_cols]
    s_true_loc = [s_ben, s_neu, s_del]

    dashlineprops = { 'lw' : SIZELINE * 2.0, 'ls' : ':', 'alpha' : 0.5, 'color' : BKCOLOR }
    histprops = dict(histtype='bar', lw=SIZELINE/2, rwidth=0.8, ls='solid', alpha=0.7, edgecolor='none',
                     orientation='horizontal')
    pprops = { 'xlim'        : [-0.04, 0.04],
               'xticks'      : [  -0.04,   -0.03,   -0.02,   -0.01,     0.,   0.01,   0.02,   0.03,   0.04],
               'yticklabels' : [r'$-4$', r'$-3$', r'$-2$', r'$-1$', r'$0$', r'$1$', r'$2$', r'$3$', r'$4$'],
               'ylim'        : [0., 0.075],
               'yticks'      : [0., 0.025, 0.05, 0.075],
               'xlabel'      : 'Frequency',
               'ylabel'      : 'Inferred selection coefficient, ' + r'$\hat{s}$' + ' (%)',
               'bins'        : np.arange(-0.04, 0.04+0.001, 0.001),
               'combine'     : True,
               'plotprops'   : histprops,
               'axoffset'    : 0.1,
               'theme'       : 'boxed' }

    for i in range(len(tags)):
        x = [np.array(df_s[cols[i]]).flatten()]
        tprops = dict(ha='left', va='center', family=FONTFAMILY, size=SIZELABEL, rotation=270, clip_on=False)
        ax_l.text(0.102, s_true_loc[i], r'$s_{%s}$' % (tags[i]), color=colors[i], **tprops)
        ax_l.axhline(y=s_true_loc[i], **dashlineprops)
        if i<len(tags)-1: mp.hist(             ax=ax_l, x=x, colors=[colors[i]], **pprops)
        else:             mp.plot(type='hist', ax=ax_l, x=x, colors=[colors[i]], **pprops)

    ## outside text labels

    tprops = dict(color=BKCOLOR, ha='center', va='center', family=FONTFAMILY, size=SIZELABEL,
                  clip_on=False, transform=fig.transFigure)
    dx = -0.04
    dy =  0.02

    ax_l.text(box_l['left']+dx,  box_l['top']+dy, 'e'.lower(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    # SAVE FIGURE

    plt.savefig('%s/revision/varying-s-%s%s' % (FIG_DIR, name, EXT), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plt.close(fig)

    print('MPL supplementary example (varying s) done.')
