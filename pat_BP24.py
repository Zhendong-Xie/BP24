import numpy as np 
from mpmath import mp
import matplotlib.pyplot as plt 
from matplotlib.font_manager import FontProperties
import pandas as pd 
from scipy.constants import N_A

import contour_plot as cplt
from constant import dsh_4T1 

adjust= [0.9, 0.15, 0.15, 0.9]
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
    h_L = 2.6  # Ligand lenth. 2.8-3 nm 

    eB = -11 # Binding energy.  <--------------------

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

    # 2. Calculate interfernce potential. 
    if h_L < 0.5 * (sigma_P / np.pi) ** 0.5:  # @ Shorter than chain distance. 
        u_j = (np.pi * h_L ** 3) / (6 * sigma_P ** (3 / 2)) * (1 - delta ** 2) ** (9 / 4)
    elif (0.5 * (sigma_P / np.pi) ** 0.5) < h_L < R:  # @ Shorter than radius. 
        u_j = ((h_L * N_PEG ** 2) * (a_PEG ** 2 / sigma_P) ** (7 / 6) / (delta ** (1 / 4) * sigma_P)) * (1 - delta ** 3)
    elif R<h_L:  # @ Longer than the radius. 
        u_j = ((R * N_PEG ** 2) * (a_PEG ** 2 / sigma_P) ** (7 / 6) / (delta ** (1 / 4) * sigma_P)) * (1 - delta ** 3)
        
    if u_j+eB>0:
        print(u_j)
    """Partition function."""
    term = 0
    # Calculate q_zeta. 
    for p in range(1, min(int(N_L), NB) + 1):
        term += mp.binomial(int(N_L),p)*mp.binomial(NB, p)*mp.factorial(p)*np.exp(-p*(eB+ u_j))
    Q_P = term
    xi = float(Q_P) 

    # Calculate the binding volume. 
    R_dm = R * 10**-8
    d_dm = h_P*10**-8 * delta
    vB = np.pi/3*((R_dm + d_dm)**3 - R_dm**3)

    # Calculate the GAG steric. 
    Ugag = np.pi *4/3 * R**3/(np.pi* dsh_4T1**2)**(3/2)
    Xi = vB*N_A* xi * np.exp(-Ugag)*np.exp(-10)

    return Xi



if __name__ == "__main__":
    fp = "./Figure/R{0}{1}.png"
    rho_P = 10**-10 # 10**6 / 0.001 l 
    print("Molarity", rho_P* N_A)
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

        # df_theta = pd.DataFrame(theta_all, columns=N_L_array)
        # pf = "./theta.xlsx"
        # df_theta.to_excel(pf, index=False)
        cticks = {"cticks":[0.00, 0.25, 0.50, 0.75, 1.00],
                  "clables":[0.00, 0.25, 0.50, 0.75, 1.00]}
        label1 = [r"$N_L$", r"$\delta_{P}$", f"$R={R_array[k]},\\theta$", None]
        cplt.contour_plot(N_L_array, delta_array, theta_all, MyColor=1, 
                          adjust=adjust, label=label1, cticks=cticks)
        # plt.ylim([0, 1])
        plt.savefig(fp.format(R_array[k],"theta"))

        label2 = [r"$N_L$", r"$\delta_{P}$", f"$R={R_array[k]},KA$", None]

        conline = [0, 24, 7]
        clim= [conline[0], conline[1]]
        cticks = {"cticks":[0, 4, 8, 12, 16, 20, 24],
                  "clables":[r"$10^{-0}$", r"$10^{-4}$", r"$10^{-8}$", 
                             r"$10^{-12}$", r"$10^{-16}$", r"$10^{-20}$",
                             r"$10^{-24}$"]}
        
        cplt.contour_plot(N_L_array, delta_array, KA_all, MyColor="Gradient",
                           label=label2, con_line=conline, clim=clim,
                            ContourLable="10^-val", adjust=adjust, cticks=cticks,
                            cbar_inv=1)
        plt.savefig(fp.format(R_array[k],"KA"))
