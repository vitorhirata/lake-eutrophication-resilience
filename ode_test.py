import random
import scipy.integrate as spi
import numpy as np

# Carpenter lake euthrophication (parameters from Dakos)
b = 0.65
r = 2.5
h = 1.95

y1=-0.1
y2=0.5
y3=1.0
y4=1.5

def possible_a(y):
    if y > y4:
        return [0]
    elif y > y3:
        return [0, 0.1]
    elif y > y2:
        return [0, 0.1, 0.2]
    elif y > y1:
        return [0, 0.1, 0.2, 0.3]

def evolve_step(y0, a, step):
    t = np.arange(0, 0 + step + 0.1, 0.1)

    sol = spi.odeint(f, y0, t, args=(a, ))
    y = sol[:, 0][-1]
    return y

def entropy(y0, a, step, time_limit, prob=1):
    possible_a_vec = possible_a(y0)
    step_prob = 1 / len(possible_a_vec)

    if time_limit == 1:
        final_prob = prob * step_prob
        return len(possible_a_vec) * ( - final_prob * np.log(final_prob))

    final_prob = prob * step_prob
    y_final = [evolve_step(y0, a, step) for a in possible_a_vec]
    return sum([entropy(y, a, step, time_limit-1, final_prob) for (a, y) in zip(possible_a_vec, y_final)])

# TODO
def entropy_iterative(y0, a, step, time_limit, prob=1):
    return 0

def f(y, t, a):
    return a - 0.65 * y + 2.5*(y**2) / (1.95**2 + y**2)

def run_scenario(t_max, y0, a, step, time_limit):
    t_final = np.arange(0, t_max + step, step)
    s_final = []
    y_final = []

    for t0 in t_final:
        y_final.append(y0)
        s_final.append(entropy(y0, a, step, time_limit))

        if t0 != 0:
            a = random.choice(possible_a(y0))

        y0 = evolve_step(y0, a, step)


    return s_final, y_final

def run_scenarios(t_max, y0, a, n_scenarious, step, time_limit):
    time_size = round((t_max + step) / step)
    s_final = np.zeros(time_size)
    y_final = np.zeros(time_size)

    for i in range(n_scenarious):
        s, y = run_scenario(t_max, y0, a, step, time_limit)
        s_final = np.sum([s_final, s], axis=0)
        y_final = np.sum([y_final, y], axis=0)

    s_final = np.divide(s_final, n_scenarious)
    y_final = np.divide(y_final, n_scenarious)

    return s_final, y_final
