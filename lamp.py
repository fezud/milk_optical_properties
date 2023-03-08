import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import pickle

P = 1
L = 0.12
d = 0.007
N = 100

D_quartz = 0.04
thickness_quartz = 0.002
thickness_milk = 0.001

milk_kappa10 = 5*10**4
quartz_kappa10 = 18.1


def unit_vector(vector):
    return vector / np.linalg.norm(vector)


def cos_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)


def circle_line_intersection(radius, point, vector, H):
    (a, b, c) = tuple(vector)
    (x0, y0, z0) = tuple(point)
    (yc, zc) = (H + d/2, 0) 

    roots = np.roots([b**2+c**2, 2*(c*(z0-zc)+b*(y0-yc)), (z0-zc)**2+(y0-yc)**2-radius**2])

    t = roots[0] if y0 + roots[0]*b < H else roots[1]

    x = x0 + t*a
    y = y0 + t*b
    z = z0 + t*c

    return (x, y, z)


def vector_length(vector):
    return np.sqrt(vector.dot(vector))

def plot_setup_zy(H, show=False):
    figure, axes = plt.subplots()
    lamp_circle = plt.Circle((0, H + d/2), radius=d/2, color='purple')
    quartz_first_layer_circle = plt.Circle((0, H + d/2), radius=D_quartz/2, color='blue', fill=False)
    quarts_second_layer_circle = plt.Circle((0, H + d/2), radius=D_quartz/2+thickness_quartz, color='blue', fill=False)
    milk_second_layer_circle = plt.Circle((0, H + d/2), radius=D_quartz/2+thickness_quartz+thickness_milk, color='red', fill=False)
    axes.add_artist(lamp_circle)
    axes.add_artist(quartz_first_layer_circle)
    axes.add_artist(quarts_second_layer_circle)
    axes.add_artist(milk_second_layer_circle)
    plt.xlim(-D_quartz/2, D_quartz/2)
    plt.ylim(0, H + d)
    plt.grid()
    if show:
        plt.show()
    return figure, axes


def plot_setup_xy(H, show=False):
    figure, axes = plt.subplots()
    axes.add_patch(Rectangle((H, -L/2), d, L, color='purple'))
    plt.plot([H + d/2 - D_quartz/2, H + d/2 - D_quartz/2], [-L/2, L/2], linestyle='-', color='blue')
    plt.plot([H + d/2 - (D_quartz/2 + thickness_quartz), H + d/2 - (D_quartz/2 + thickness_quartz)], [-L/2, L/2], linestyle='-', color='blue')
    plt.plot([H + d/2 - (D_quartz/2 + thickness_quartz + thickness_milk), H + d/2 - (D_quartz/2 + thickness_quartz + thickness_milk)], [-L/2, L/2], linestyle='-', color='red')
    plt.xlim(0, H + d + d/5)
    plt.ylim(-(L/2 + d/10), L/2 + d/10)
    plt.grid()
    if show:
        plt.show()
    return figure, axes   



def Ps(H, initial_x=0):
    phi_max = np.arccos((d/2)/((d/2) + H))
    # максимальный угол, излучение с которого
    # будет попадать в точку (0, 0, 0)
    d_phi = phi_max/N  # шаг изменения угла
    dx = L/N  # шаг изменения координаты вдоль оси лампы

    zero_vector = np.array([0, 1, 0])  # нормальный вектор в точке (0, 0, 0)

    dPs = P/(np.pi*d*L)
    dP = dx*d_phi*(d/2)*dPs/(np.pi)  # излучение с элементарной площадки

    Ps_zero = 0  # суммарное излучение, полученное точкой

    for x in np.arange(-L/2, L/2 + dx, dx):
        for phi in np.arange(0, phi_max + d_phi, d_phi):
            y_r = np.cos(phi)*(d/2)
            z_r = np.sin(phi)*(d/2)

            y = H + d/2 - y_r
            z = z_r

            length_vector = np.array([x, y, z]) - np.array([initial_x, 0, 0])


            (x_quartz_intersection, y_quartz_intersection, z_quartz_intersection) = circle_line_intersection(D_quartz/2, [x, y, z], length_vector, H)
            quartz_intersection_vector = np.array([x, y, z]) - np.array([x_quartz_intersection, y_quartz_intersection, z_quartz_intersection])

            # plt.plot([z_quartz_intersection, z], [y_quartz_intersection, y], 'bo', linestyle="--", color='yellow')
            

            (x_milk_intersection, y_milk_intersection, z_milk_intersection) = circle_line_intersection(D_quartz/2 + thickness_quartz, [x, y, z], length_vector, H)
            milk_intersection_vector = np.array([x, y, z]) - np.array([x_milk_intersection, y_milk_intersection, z_milk_intersection])
            

            # вектор, соединяющий точку на лампе с координатами (x, y, z)
            # с точкой (0, 0, 0) (вектор "длины")
            radius_vector = np.array([x, y, z]) - np.array([x, H + d/2, 0])
            

            # вектор, соединяющий точку на лампе с координатами (x, y, z)
            # с точкой на оси лампы с координатами (x, H + d/2, 0)

            (length, quartz_length, milk_length) = (vector_length(vector) for vector in [length_vector, quartz_intersection_vector, milk_intersection_vector])
            # расстояние до начальной точки P.S. .dot - скалярное умножение

            zero_length_cos = cos_between(length_vector, zero_vector)
            # косинус между нормалью в начальной точке и вектором "длины"
            radius_length_cos = cos_between(-length_vector, radius_vector)
            # косинус между нормалью точке на лампе и вектором "длины"

            deltaR_milk = length - milk_length
            deltaR_quartz = length - (deltaR_milk + quartz_length)

            cos_multiplied = zero_length_cos*radius_length_cos

            milk_absorption = 10 ** (-milk_kappa10*deltaR_milk)
            quartz_absorption = 10 ** (-quartz_kappa10*deltaR_quartz)

            Ps_zero += (dP/(length*length)) * \
                cos_multiplied * 2 * \
                milk_absorption * \
                quartz_absorption
            # 2 – из симметрии

            # if x >= 10*dx and x <= L/2:
            #     plot_setup_zy(H)
            #     plt.plot([0, z], [0, y], 'bo', linestyle="--", color='black')
            #     plt.plot([z_quartz_intersection], [y_quartz_intersection], marker="o", color='green')
            #     plt.plot([z_milk_intersection], [y_milk_intersection], marker="o", color='green')
            #     plt.plot([z_milk_intersection, z], [y_milk_intersection, y], 'bo', linestyle="--", color='yellow')
            #     plt.plot([0, z], [H + d/2, y], 'bo', linestyle="--", color='black')
            #     plt.show()

            # if x >= 10*dx and x <= L/2:
            #     plot_setup_xy(H)
            #     plt.plot([0, y], [0, x], 'bo', linestyle="--", color='black')
            #     plt.plot([y_quartz_intersection], [x_quartz_intersection], marker="o", color='green')
            #     plt.plot([y_milk_intersection], [x_milk_intersection], marker="o", color='green')
            #     plt.plot([y_milk_intersection, y], [x_milk_intersection, x], 'bo', linestyle="--", color='yellow')
            #     plt.plot([H + d/2, y], [x, x], 'bo', linestyle="--", color='black')
            #     plt.show()
    return Ps_zero


# H = np.arange(D_quartz/2 + thickness_quartz - d/2, D_quartz/2 + thickness_quartz + thickness_milk - d/2, thickness_milk/20)
H = np.arange(0.01850, 0.01860, 0.000001)
fractions = [Ps(h, initial_x=0) for h in H]
# for h in H:
#     fractions.append(Ps(h, initial_x=0))


# # print(Ps(D_quartz/2 + thickness_quartz - d/2, initial_x=0))

# print(pickle.dumps(H))
# with open('data.pickle', 'wb') as file:
#     data = {'H, m': H, 'P, отн. ед': fractions}
#     pickle.dump(data, file)

with open('data.pickle', 'rb') as file:
    h: np.arange = pickle.load(file)
    print(h)

# plt.figure(figsize=(15, 13), facecolor='white')
# plt.plot(H[:], fractions[:], color='tab:blue')
# plt.show()
# plt.grid()
# plt.yscale("log")
# plt.xlabel('H, m')
# plt.ylabel('P, отн. ед.')
