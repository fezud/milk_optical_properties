import pickle
import numpy as np
import matplotlib.pyplot as plt

D90 = 0.22 * (10**-3) * (10**4)  # изначально 0.22 мДж/см2


with open('data.pickle', 'rb') as file:
    data: dict = pickle.load(file)
    H: np.arange = data['H, m']
    fractions: list = data['P, отн. ед']

    h = [(value-H[0])*10**6 for value in H]

    eta = [(1 - 10**(-value*2*3/D90)) for value in fractions]

    H_full_layer = np.arange(0.01850, 0.0195, 0.000001)
    H_full_layer = [(h - H_full_layer[0])*10**6 for h in H_full_layer]

    eta_full_layer = np.zeros(len(H_full_layer))
    eta_full_layer[0:len(eta)] = eta
    eta_full_layer[len(eta):] = [0 for _ in range(len(eta), len(eta_full_layer))]
    print(eta_full_layer)

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

    plt.figure(figsize=(15, 13), facecolor='white')
    plt.plot(H_full_layer[:], eta_full_layer[:], color='tab:red')
    plt.xlabel('H_full_layer, мкм')
    plt.ylabel('степень витаминизации, 1')
    plt.grid()
    # plt.show()

    numerator = 0
    denominator = 0

    for i in range(2, len(eta_full_layer)):
        numerator += eta_full_layer[i]*(H_full_layer[i]+0.01850)*(H_full_layer[i]-H_full_layer[i-1])
        denominator += (H_full_layer[i]+0.01850)*(H_full_layer[i]-H_full_layer[i-1])

    average_eta = numerator/denominator
    print(f'Eta average value = {average_eta}')

    n = np.log(0.5)/np.log(1-average_eta)
    print(n)


