import pickle
import numpy as np
import matplotlib.pyplot as plt

D90 = 0.22 * (10**-3) * (10**4) # изначально 0.22 мДж/см2


with open('data.pickle', 'rb') as file:
    data: dict = pickle.load(file)
    H: np.arange = data['H, m']
    fractions: list = data['P, отн. ед']

    h = [(value-H[0])*10**6 for value in H]

    eta = [(1 - 10**(-value*12/D90)) for value in fractions]

    plt.figure(figsize=(15, 13), facecolor='white')
    plt.plot(h[:], fractions[:], color='tab:blue')
    plt.xlabel('H, мкм')
    plt.ylabel('P, отн. ед.')
    plt.grid()

    plt.figure(figsize=(15, 13), facecolor='white')
    plt.plot(h[:], eta[:], color='tab:blue')
    plt.xlabel('H, мкм')
    plt.ylabel('степень витаминизации, 1')
    plt.grid()
    plt.show()