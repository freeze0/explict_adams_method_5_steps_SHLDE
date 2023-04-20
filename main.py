import numpy as np
import csv

def adams5(f, tspan, y0, h):
    # f - функция, задающая правую часть уравнений
    # tspan - временной интервал для решения (предполагается, что tspan[0] - начальный момент времени,
    # а tspan[1] - конечный момент времени)
    # y0 - начальные значения вектора решений
    # h - шаг интегрирования

    # Определяем количество шагов интегрирования
    n = int(np.round((tspan[1] - tspan[0]) / h))

    # Создаем массивы для хранения решений
    t = np.linspace(tspan[0], tspan[1], n + 1, dtype=np.float16)
    y = np.zeros((n + 1, len(y0)), dtype=np.float16)

    # Заполняем первые 5 точек используя метод Рунге-Кутта 4 порядка
    y[0, :] = y0
    for i in range(4):
        k1 = f(t[i], y[i, :])
        k2 = f(t[i] + h / 2, y[i, :] + h / 2 * k1)
        k3 = f(t[i] + h / 2, y[i, :] + h / 2 * k2)
        k4 = f(t[i] + h, y[i, :] + h * k3)
        y[i + 1, :] = y[i, :] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    # Используем пятишаговый метод Адамса для остальных точек
    for i in range(4, n):
        fn = f(t[i], y[i, :])
        fn_1 = f(t[i - 1], y[i - 1, :])
        fn_2 = f(t[i - 2], y[i - 2, :])
        fn_3 = f(t[i - 3], y[i - 3, :])
        fn_4 = f(t[i - 4], y[i - 4, :])
        y[i + 1, :] = y[i, :] + h / 720 * (1901 * fn - 2774 * fn_1 + 2616 * fn_2 - 1274 * fn_3 + 251 * fn_4)

    return t, y


# # Задаем функцию, задающую правую часть системы ZADACHA 1
# def f(t, y):
#     dx_dt = np.zeros(y.shape)
#     dx_dt[0] = 2*y[0]-y[1]-y[2]
#     dx_dt[1] = 3*y[0]-2*y[1]-3*y[2]
#     dx_dt[2] = -y[0]+y[1]+2*y[2]
#     return dx_dt
#
# def exact(t):
#     y1 = 2*np.exp(t)
#     y2 = np.exp(t)
#     y3 = np.exp(t)
#     return [y1, y2, y3]
#
#
# # Начальные значения и временной интервал для решения
# y0 = np.array([2, 1, 1])
# tspan = [0, 3]
#
# # Решение системы с помощью пятишагового метода Адамса с шагом 0.1
# x_h, y_adams_h = adams5(f, tspan, y0, 0.1)
# y_h = exact(x_h)
# y_h = np.reshape(y_h, (3, 31))
#
# # Решение системы с помощью пятишагового метода Адамса с шагом 0.05
# x_h2, y_adams_h2 = adams5(f, tspan, y0, 0.05)
# y_h2 = exact(x_h2)
# y_h2 = np.reshape(y_h2, (3, 61))
#
# with open('answer.csv', 'w', newline='') as csvfile:
#     fieldnames = ['x_h', 'y1_adams_h', 'y2_adams_h', 'y3_adams_h', 'y1_toch_h', 'y2_toch_h', 'y3_toch_h', 'delta1_h',
#                   'delta2_h', 'delta3_h', 'x_h2', 'y1_adams_h2', 'y2_adams_h2', 'y3_adams_h2', 'y1_toch_h2',
#                   'y2_toch_h2', 'y3_toch_h2', 'delta1_h2', 'delta2_h2', 'delta3_h2']
#     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
#     writer.writeheader()
#     for i in range(len(x_h2)):
#         writer.writerow({'x_h2': x_h2[i], 'y1_adams_h2': y_adams_h2[i][0], 'y2_adams_h2': y_adams_h2[i][1],
#                         'y3_adams_h2': y_adams_h2[i][2], 'y1_toch_h2': y_h2[0][i], 'y2_toch_h2': y_h2[1][i],
#                          'y3_toch_h2': y_h2[2][i], 'delta1_h2': abs(y_adams_h2[i][0] - y_h2[0][i]),
#                          'delta2_h2': abs(y_adams_h2[i][1] - y_h2[1][i]), 'delta3_h2': abs(y_adams_h2[i][2] - y_h2[2][i])})
#     for i in range(len(x_h)):
#         writer.writerow({'x_h': x_h[i], 'y1_adams_h': y_adams_h[i][0], 'y2_adams_h': y_adams_h[i][1],
#                          'y3_adams_h': y_adams_h[i][2], 'y1_toch_h': y_h[0][i], 'y2_toch_h': y_h[1][i],
#                          'y3_toch_h': y_h[2][i], 'delta1_h': abs(y_adams_h[i][0] - y_h[0][i]),
#                          'delta2_h': abs(y_adams_h[i][1] - y_h[1][i]), 'delta3_h': abs(y_adams_h[i][2] - y_h[2][i])})
# Задаем функцию, задающую правую часть системы

def f(t, y):
    dx_dt = np.zeros(y.shape)
    dx_dt[0] = y[1]**2-2*t*y[1]-2*y[1]-y[0]
    dx_dt[1] = 2*y[0]+2*t**2+np.exp(2*t-2*y[1])
    return dx_dt

def exact(t):
    y1 = -(t**2)
    y2 = t
    return [y1, y2]


# Начальные значения и временной интервал для решения
y0 = np.array([0, 0])
tspan = [0, 3.5]

# Решение системы с помощью пятишагового метода Адамса с шагом 0.1
x_h, y_adams_h = adams5(f, tspan, y0, 0.1)
y_h = exact(x_h)
y_h = np.reshape(y_h, (2, len(y_adams_h)))

# Решение системы с помощью пятишагового метода Адамса с шагом 0.05
x_h2, y_adams_h2 = adams5(f, tspan, y0, 0.05)
y_h2 = exact(x_h2)
y_h2 = np.reshape(y_h2, (2, len(y_adams_h2)))

with open('answer.csv', 'w', newline='') as csvfile:
    fieldnames = ['x_h', 'y1_adams_h', 'y2_adams_h', 'y1_toch_h', 'y2_toch_h', 'delta1_h', 'delta2_h', 'x_h2',
                  'y1_adams_h2', 'y2_adams_h2', 'y1_toch_h2', 'y2_toch_h2', 'delta1_h2', 'delta2_h2']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for i in range(len(x_h2)):
        writer.writerow({'x_h2': x_h2[i], 'y1_adams_h2': y_adams_h2[i][0], 'y2_adams_h2': y_adams_h2[i][1],
                         'y1_toch_h2': y_h2[0][i], 'y2_toch_h2': y_h2[1][i],
                         'delta1_h2': abs(y_adams_h2[i][0] - y_h2[0][i]),
                         'delta2_h2': abs(y_adams_h2[i][1] - y_h2[1][i])})
    for i in range(len(x_h)):
        writer.writerow({'x_h': x_h[i], 'y1_adams_h': y_adams_h[i][0], 'y2_adams_h': y_adams_h[i][1],
                         'y1_toch_h': y_h[0][i], 'y2_toch_h': y_h[1][i], 'delta1_h': abs(y_adams_h[i][0] - y_h[0][i]),
                         'delta2_h': abs(y_adams_h[i][1] - y_h[1][i])})
