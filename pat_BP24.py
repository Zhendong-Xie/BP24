import numpy as np 
from mpmath import mp
import matplotlib.pyplot as plt 
from matplotlib.font_manager import FontProperties
import pandas as pd 
from scipy.constants import N_A

import contour_plot as cplt

adjust = [0.9,0.15,0.17,0.95]
# Here you should give the polymerization, monomer size, distance of the PEG.
def pat_bp24(N_PEG, R, delta, N_L, NB=3, output=None):
    """
    This function is calculate the partition function between the NP and amyloid beta.

    Parameters:
    -----------
    R: radius. (nm).

    delta: interference parameter, dimensionless.

    N_L: number of ligand per particle.
    """

    # 1. Particle properties.
    a_PEG = 0.35 # PEG size (radius).
    d_PEG = 2.5 # Distance between 2 PEG <-------------------- modify parmeter.
    h_L = 2.6 + 0.6  # Ligand lenth.
    eB = -11.6 # Binding energy.  <--------------------

    # 2. Calculate the insertion steric.
    """Don't change the compress steric."""
    # 2.1 PEG brush height on flat surface.
    h_0 = N_PEG * (( 4*(a_PEG**3) * a_PEG ** 2) / (9* (d_PEG ** 2))) ** (1 / 3)
    # Geometry parameter gamma.
    gamma = (h_0 / R + 1) ** 2 if h_0 / R <= (np.sqrt(3) - 1) else 3
    # Height on curve surface.
    h_P = R * ((1 + ((gamma + 2) * h_0) / (3 * R)) ** (3 / (gamma + 2)) - 1)
    # Density on the curve surface.
    sigma_P = np.pi * d_PEG ** 2 * (1 + (delta * h_P / R)) ** (gamma - 1)

    # 2.2. Calculate interfernce potential.
    # @ Shorter than chain distance.
    if h_L < 0.5 * (sigma_P / np.pi) ** 0.5:
        u_j=(np.pi*h_L**3)/(6*sigma_P**(3/2))*(1-delta**2)**(9/4)
    # Shorter than radius.
    elif (0.5 * (sigma_P / np.pi) ** 0.5) < h_L < R:
        denominator = (delta ** (1 / 4) * sigma_P)
        numerator = (h_L * N_PEG ** 2) * (a_PEG ** 2 / sigma_P) ** (7 / 6)
        u_j = numerator/ denominator * (1 - delta ** 3)
    # @ Longer than the radius.
    elif R<h_L:
        denominator = (delta ** (1 / 4) * sigma_P)
        numerator = (R * N_PEG ** 2) * (a_PEG ** 2 / sigma_P) ** (7 / 6)
        u_j = numerator/ denominator * (1 - delta ** 3)

    """Partition function."""
    term = 0
    u_s = 0
    # Calculate q_zeta.
    for p in range(1, min(int(N_L), NB) + 1):
        # u_s = (int(N_L) - p) * A *  delta ** 3   # The interference steric from unbinding ligand.
        term += mp.binomial(int(N_L),p)*mp.binomial(NB, p)*mp.factorial(p)*np.exp(-p*(eB+ u_j))*np.exp(-u_s)
    Q_P = term
    xi = float(Q_P)

    # Calculate the binding volume.
    R_dm = R * 10**-8
    d_dm = h_P*10**-8 * delta
    vB = np.pi/3*((R_dm + d_dm)**3 - R_dm**3)

    # Calculate the GAG steric.
    Ugag = 5
    Xi = vB*N_A* xi * np.exp(-Ugag)

    return Xi



if __name__ == "__main__":
    fp = "./Figure/R{0}{1}.png"
    rho_P = 10**-10 # 10**6 / 0.001 l
    print("Molarity", rho_P/ N_A)
    R_array = [6.5, 10.5, 19.5]
    N_PEG_array = [46, 114, 228]

    #                       NPEG, R, Delta
    para_array = np.array([[46, 6.5, 0.58],
                           [46, 6.5, 0.82],
                           [114, 10.5, 0.28],
                           [114, 10.5, 0.39],
                           [114, 10.5, 0.66],
                           [114, 10.5, 1],
                           [228, 19.5, 0.16],
                           [228, 19.5, 0.23],
                           [228, 19.5, 0.38],
                           [228, 19.5, 0.58],
                           [228, 19.5, 1],
                           ])
    # The formulation of NL.
    N_L_exp = np.array([5, 10, 20, 40])

    delta_array = np.linspace(0.01, 1, 100)
    N_L_array = np.linspace(1, 40, 40)

    row = len(delta_array)
    col = len(N_L_array)

    for k in range(len(R_array)):
        theta_all = np.zeros([row, col])
        KA_all = np.zeros([row, col])
        """The binding fraction versus the delta and NL."""
        for x_i in range(row):
            for x_j in range(col):
                xi_i = pat_bp24(N_PEG_array[k], R_array[k], delta_array[x_i], N_L_array[x_j])
                theta_i = rho_P *xi_i / (1 +  rho_P * xi_i )
                theta_all[x_i, x_j] = theta_i

                KA_all[x_i, x_j] = np.log10(xi_i)
        print(KA_all)
        # df_theta = pd.DataFrame(theta_all, columns=N_L_array)
        # pf = "./theta.xlsx"
        # df_theta.to_excel(pf, index=False)
        # theta plot
        label1 = [r"$N_L$", r"$\delta_{P}$", f"$R={R_array[k]},\\theta$", None]
        cplt.contour_plot(N_L_array, delta_array, theta_all, MyColor=1, label=label1,adjust=adjust)
        # plt.ylim([0, 1])
        plt.savefig(fp.format(R_array[k],"theta"))

        label2 = [r"$N_L$", r"$\delta_{P}$", f"$R={R_array[k]},KA$", None]
        # avidity plot
        conline = [0, 24, 7]
        clim= [conline[0], conline[1]]
        cplt.contour_plot(N_L_array, delta_array, KA_all, MyColor="Gradient",
                           label=label2, con_line=conline, clim=clim,
                            ContourLable="10^-val",adjust=adjust)
        plt.savefig(fp.format(R_array[k],"KA"))
