# -*- coding: utf-8 -*-
"""

Code to generate Fig. 2D.

Fig. 2D shows the approximate and exact loss ratio between the SI and MI motifs
defined in the analytical study of information loss.
The three "idealized" motifs are named "SI(A)", "SI(B)" and "MI", and their
information losses are studied as a function of N (the number of sites in each motif).

"""

from math import log2
import pandas as pd
import matplotlib.pyplot as plt


def H(probs):
    ''' Returns the Shannon Entropy of the probability distribution `probs`.'''
    return -sum(p * log2(p) for p in probs if p > 0)


def expected_loss_SI_A(N):
    '''
    Returns the expected information loss due to a point mutation for Motif "SI(A)".
    `N` is the number of sites.
    '''
    # One point mutation changes (two) frequencies by +/- (1/N)
    eps = 1 / N

    # Case 1: mutate X
    # X: A → other letter
    Hx_before = 0  # because fully conserved
    # After mutation
    Hx_after = H([1 - eps, eps, 0, 0])
    loss_X = Hx_after - Hx_before
    
    # Case 2: mutate Y
    # Y was uniform; one mutation changes one count by 1
    Hy_before = 2  # because uniform
    # After mutation
    Hy_after = H([0.25 + eps, 0.25 - eps, 0.25, 0.25])
    loss_Y = Hy_after - Hy_before
    
    return 0.5 * (loss_Y + loss_X)


def expected_loss_SI_B(N):
    '''
    Returns the expected information loss due to a point mutation for Motif "SI(B)".
    `N` is the number of sites.
    '''
    # One point mutation changes (two) frequencies by +/- (1/N)
    eps = 1 / N
    
    # Case 1: mutate X
    Hx_before = 1  # because half conserved
    
    # Mutation stays within {A,C}
    Hx_after = H([0.5 - eps, 0.5 + eps, 0, 0])
    loss_X_toC = Hx_after - Hx_before
    
    # Mutation introduces new letter
    Hx_after = H([0.5 - eps, 0.5, eps, 0])
    loss_X_toGT = Hx_after - Hx_before
    
    loss_X = ((1/3)*loss_X_toC) + ((2/3)*loss_X_toGT)
    
    # Case 2: mutate Y
    # Same by symmetry --> avg loss = avg loss_X
    return loss_X


def expected_loss_MI(N):
    '''
    Returns the expected information loss due to a point mutation for Motif "MI".
    `N` is the number of sites.
    '''
    # One point mutation changes (two) frequencies by +/- (1/N)
    eps = 1 / N
    
    # Calculate MI as MIxy = Hx + Hy - Hxy
    MI_before = 2
    
    # Case 1: mutate X
    Hx = H([0.25+eps, 0.25-eps, 0.25, 0.25])
    Hy = 2  # because uniform
    # Four dimers are now almost equally likely; and new dimer has frequency eps
    Hxy = H([0.25-eps, 0.25, 0.25, 0.25, eps] + [0]*11)
    MI_after = Hx + Hy - Hxy
    
    loss = MI_before - MI_after
    
    # Case 2: mutate Y
    # Same by symmetry --> avg loss = loss
    # So we can just return the loss from Case 1
    return loss


def compute_information_losses(exp_min=2, exp_max=50):
    '''
    Computes the approximate (analytically derived) and exact information loss
    for motifs "SI(A)", "SI(B)" and "MI", as well as the loss ratio between
    SI(A) and MI, and between SI(B) and MI.
    
    The calculation is performed for all values of N (= number of sites) from
    2^exp_min to 2^exp_max.
    '''
    
    exponents = list(range(2, 50))
    Nvals = [2**e for e in exponents]
    
    # Analytical approximation of losses
    # ----------------------------------
    
    SI_A_approx = [log2(N)/(2*N) for N in Nvals]
    SI_B_approx = [2 * log2(N)/(3*N) for N in Nvals]
    MI_approx = [log2(N)/N for N in Nvals]
    
    # Compute exact losses
    # --------------------
    
    # Losses
    loss_SI_A = []
    loss_SI_B = []
    loss_MI = []
    
    # Loss Ratios
    r_SIA_MI = []
    r_SIB_MI = []
    
    for exp in exponents:
        N = 2**exp
        
        # Information losses
        l_siA = expected_loss_SI_A(N)
        l_siB = expected_loss_SI_B(N)
        l_mi = expected_loss_MI(N)
        
        loss_SI_A.append(l_siA)
        loss_SI_B.append(l_siB)
        loss_MI.append(l_mi)
        
        # Ratios
        ratio_SIA_MI = l_siA / l_mi
        ratio_SIB_MI = l_siB / l_mi
        
        r_SIA_MI.append(ratio_SIA_MI)
        r_SIB_MI.append(ratio_SIB_MI)
        
        # print(f"N = {N}")
        # print(f"  SI_A loss = {l_siA:.6e}")
        # print(f"  SI_B loss = {l_siB:.6e}")
        # print(f"  MI   loss = {l_mi:.6e}")
        # print(f"  SI_A / MI ≈ {ratio_SIA_MI:.3f}")
        # print(f"  SI_B / MI ≈ {ratio_SIB_MI:.3f}")
        # print()
    
    data = pd.DataFrame({'N': Nvals,
                         'SI_A_approx': SI_A_approx,
                         'SI_B_approx': SI_B_approx,
                         'MI_approx': MI_approx,
                         'loss_SI_A': loss_SI_A,
                         'loss_SI_B': loss_SI_B,
                         'loss_MI': loss_MI,
                         'r_SIA_MI': r_SIA_MI,
                         'r_SIB_MI': r_SIB_MI})
    return data


def plot_losses(data, outfilename=None):
    # Plot LOSSES
    RIGHT_BOUND = 5
    Nvals_small = data['N'][:RIGHT_BOUND]
    plt.plot(Nvals_small, data['loss_SI_A'][:RIGHT_BOUND], color='blue', label=r'$\mathbb{E}[\Delta_{SI(A)}]$')
    plt.plot(Nvals_small, data['SI_A_approx'][:RIGHT_BOUND], linestyle='dashed', color='blue', label=r'$\mathregular{r_{SI(A)}}$')
    plt.plot(Nvals_small, data['loss_SI_B'][:RIGHT_BOUND], color='orange', label=r'$\mathbb{E}[\Delta_{SI(B)}]$')
    plt.plot(Nvals_small, data['SI_B_approx'][:RIGHT_BOUND], linestyle='dashed', color='orange', label=r'$\mathregular{r_{SI(B)}}$')
    plt.plot(Nvals_small, data['loss_MI'][:RIGHT_BOUND], color='green', label=r'$\mathbb{E}[\Delta_{MI}]$')
    plt.plot(Nvals_small, data['MI_approx'][:RIGHT_BOUND], linestyle='dashed', color='green', label=r'$\mathregular{r_{MI}}$')
    plt.xscale("log", base=2)
    plt.xlabel('N')
    plt.ylabel('information loss (bits)')
    plt.grid()
    plt.legend()
    if outfilename is not None:
        plt.savefig(outfilename, bbox_inches='tight')
    plt.show()


def plot_losses_convergence(data, outfilename=None):
    # Plot convergence of losses
    Nvals = data['N']
    conv_SI_A = [data['loss_SI_A'][i]/data['SI_A_approx'][i] for i in range(len(Nvals))]
    conv_SI_B = [data['loss_SI_B'][i]/data['SI_B_approx'][i] for i in range(len(Nvals))]
    conv_MI   = [data['loss_MI'][i]/data['MI_approx'][i] for i in range(len(data['N']))]
    plt.plot(Nvals, conv_SI_A, color='blue', label=r'$\mathbb{E}[\Delta_{SI(A)}]/r_{SI(A)}$')
    plt.plot(Nvals, conv_SI_B, color='orange', label=r'$\mathbb{E}[\Delta_{SI(B)}]/r_{SI(B)}$')
    plt.plot(Nvals, conv_MI, color='green', label=r'$\mathbb{E}[\Delta_{MI}]/r_{MI}$')
    plt.hlines(1, min(Nvals), max(Nvals), color='black', linestyles='dashed')
    plt.xscale("log", base=2)
    plt.xlabel('N')
    plt.grid()
    plt.legend()
    if outfilename is not None:
        plt.savefig(outfilename, bbox_inches='tight')
    plt.show()


def plot_loss_ratios(data, outfilename=None, ax=None):
    # Plot LOSS RATIOS
    Nvals = data['N']
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(Nvals, data['r_SIA_MI'], color='blue', label=r'$\mathbb{E}[\Delta_{SI(A)}]$ / $\mathbb{E}[\Delta_{MI}]$')
    ax.hlines(1/2, min(Nvals), max(Nvals), color='blue', linestyles='dashed', label='$r_{SI(A)}/r_{MI}=1/2$')
    ax.plot(Nvals, data['r_SIB_MI'], color='orange', label=r'$\mathbb{E}[\Delta_{SI(B)}]$ / $\mathbb{E}[\Delta_{MI}]$')
    ax.hlines(2/3, min(Nvals), max(Nvals), color='orange', linestyles='dashed', label='$r_{SI(B)}/r_{MI}=2/3$')
    ax.set_ylim((0, 1))
    ax.set_xlim((min(Nvals), max(Nvals)))
    ax.set_xscale("log", base=2)
    ax.set_xlabel('N')
    ax.set_ylabel('SI loss / MI loss')
    ax.grid(True)
    ax.legend(loc='lower right')
    if outfilename is not None:
        ax.get_figure().savefig(outfilename, bbox_inches='tight')
    return ax


def plot_loss_ratios_convergence(data, outfilename=None):
    # Plot convergence of loss ratios
    Nvals = data['N']
    r_siA_mi_conv = [r/(1/2) for r in data['r_SIA_MI']]
    r_siB_mi_conv = [r/(2/3) for r in data['r_SIB_MI']]
    plt.plot(Nvals, r_siA_mi_conv, color='blue', label=r'$\mathregular{r_{SI(A)}/(1/2)}$')
    plt.plot(Nvals, r_siB_mi_conv, color='orange', label=r'$\mathregular{r_{SI(B)}/(2/3)}$')
    plt.hlines(1, min(Nvals), max(Nvals), color='black', linestyles='dashed')
    plt.ylim((0.5, 1.5))
    plt.xscale("log", base=2)
    plt.xlabel('N')
    plt.grid()
    plt.legend()
    if outfilename is not None:
        plt.savefig(outfilename, bbox_inches='tight')
    plt.show()

def make_figure_2D(data=None, ax=None, outfilename='Fig_2D.png'):
    ''' Generates Figure 2D (loss ratios plot). '''
    if data is None:
        data = compute_information_losses()
    plot_loss_ratios(data, outfilename=outfilename, ax=ax)
    return data


def main():
    data = compute_information_losses()
    plot_losses(data)
    plot_losses_convergence(data)
    plot_loss_ratios(data)
    plot_loss_ratios_convergence(data)
    make_figure_2D(data)


if __name__ == "__main__":
    main()














