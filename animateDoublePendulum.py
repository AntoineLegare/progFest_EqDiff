import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


def d2theta_1(theta_1, theta_2, dtheta_1, dtheta_2, params):
    g, L_1, L_2, m_1, m_2 = params
    term1 = -g * (2 * (m_1 + m_2)) * np.sin(theta_1)
    term2 = -m_2 * g * np.sin(theta_1 - (2 * theta_2))
    term3 = -2 * np.sin(theta_1 - theta_2) * m_2 * (((dtheta_2 ** 2) * L_2) + ((dtheta_1 ** 2) * L_1 * np.cos(theta_1 - theta_2)))
    denominator = L_1 * ((2 * m_1) + m_2 - (m_2 * np.cos((2 * theta_1) - (2 * theta_2))))
    return (term1 + term2 + term3) / denominator


def d2theta_2(theta_1, theta_2, dtheta_1, dtheta_2, params):
    g, L_1, L_2, m_1, m_2 = params
    factor = 2 * np.sin(theta_1 - theta_2)
    term1 = (dtheta_1 ** 2) * L_1 * (m_1 + m_2)
    term2 = g * (m_1 + m_2) * np.cos(theta_1)
    term3 = (dtheta_2 ** 2) * L_2 * m_2 * np.cos(theta_1 - theta_2)
    denominator = L_2 * ((2 * m_1) + m_2 - (m_2 * np.cos((2 * theta_1) - (2 * theta_2))))
    return (factor * (term1 + term2 + term3)) / denominator


def solveDoublePendulum(duration, step, initialConditions, params):
    theta_1_0, theta_2_0, dtheta_1_0, dtheta_2_0 = initialConditions

    N = int(duration / step)
    t = np.linspace(0, duration, N + 1)

    theta_1, theta_2 = [theta_1_0], [theta_2_0]
    dtheta_1, dtheta_2 = [dtheta_1_0], [dtheta_2_0]

    for i in range(N):
        dtheta_1.append(dtheta_1[-1] + step * d2theta_1(theta_1[-1], theta_2[-1], dtheta_1[-1], dtheta_2[-1], params))
        dtheta_2.append(dtheta_2[-1] + step * d2theta_2(theta_1[-1], theta_2[-1], dtheta_1[-1], dtheta_2[-1], params))
        theta_1.append(theta_1[-1] + step * dtheta_1[-1])
        theta_2.append(theta_2[-1] + step * dtheta_2[-1])

    return t, theta_1, theta_2


def getCartesianCoordinates(theta_1_list, theta_2_list, L_1, L_2):
    x_1, y_1, x_2, y_2 = [], [], [], []
    for i, theta_1 in enumerate(theta_1_list):
        theta_2 = theta_2_list[i]
        x_1.append(-L_1 * np.sin(theta_1))
        y_1.append(-L_1 * np.cos(theta_1))
        x_2.append(-L_2 * np.sin(theta_2) + x_1[-1])
        y_2.append(-L_2 * np.cos(theta_2) + y_1[-1])
    return x_1, y_1, x_2, y_2


def animateDoublePendulum(theta_1, theta_2, L_1, L_2, step):
    x_1, y_1, x_2, y_2 = getCartesianCoordinates(theta_1, theta_2, L_1, L_2)

    fig, ax = plt.subplots(figsize=(5, 5))
    line1, = ax.plot([0, x_1[0]], [0, y_1[0]], linewidth=2, color='black', zorder=10)
    mass1 = ax.scatter(x_1[0], y_1[0], s=75, linewidth=2, color='black', zorder=10)
    line2, = ax.plot([x_1[0], x_2[0]], [y_1[0], y_2[0]], linewidth=2, color='black', zorder=10)
    mass2 = ax.scatter([x_2[0]], [y_2[0]], s=75, linewidth=2, color='black', zorder=10)
    trajectory, = ax.plot(x_2[0], y_2[0], linewidth=1, color='red', zorder=0)
    ax.set_xlim([-1.1 * (L_1 + L_2), 1.1 * (L_1 + L_2)])
    ax.set_ylim([-1.1 * (L_1 + L_2), 1.1 * (L_1 + L_2)])
    plt.axis('off')


    frameStep = int(0.01 / step)
    def drawframe(i):
        i *= frameStep
        line1.set_data([0, x_1[i]], [0, y_1[i]])
        mass1.set_offsets(np.array([x_1[i], y_1[i]]))
        line2.set_data([x_1[i], x_2[i]], [y_1[i], y_2[i]])
        mass2.set_offsets(np.array([x_2[i], y_2[i]]))
        trajectory.set_data(x_2[:i], y_2[:i])
        return line1, mass1, line2, mass2, trajectory,

    anim = animation.FuncAnimation(fig, drawframe, frames=int(len(x_1)/frameStep), interval=10, blit=True)
    plt.show()


if __name__ == '__main__':

    g = 9.8  # m/s^2
    L_1 = 3  # m
    L_2 = 1  # m
    m_1 = 2  # kg
    m_2 = 1  # kg
    params = [g, L_1, L_2, m_1, m_2]

    theta_1_0 = 4
    theta_2_0 = 3
    dtheta_1_0 = 0
    dtheta_2_0 = 2
    initialConditions = [theta_1_0, theta_2_0, dtheta_1_0, dtheta_2_0]

    duration = 100
    step = 0.01
    t, theta_1, theta_2 = solveDoublePendulum(duration, step, initialConditions, params)

    animateDoublePendulum(theta_1, theta_2, L_1, L_2, step)