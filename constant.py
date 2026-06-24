import numpy as np 

dsh_4T1 = 15 

# CD24 per cell 
N_CD24_4T1 = 3.82 * 10**5 

# surface of 4T1 cell
S = 4*np.pi * 12**2  
sigmma_4T1 = N_CD24_4T1 / S  # N_CD24 per um2 


if __name__ =="__main__":
    print("Surface area of 4T1", S)
    print("Density of CD24 um2", np.sqrt(sigmma_4T1))
    R = 0.019
    NB_19 = np.pi * R**2 * sigmma_4T1 
    print(NB_19) 