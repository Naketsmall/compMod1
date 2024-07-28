import numpy as np


def Wr(field, n, i, tau, h, f_a):
    return (f_a(n * tau, i * h) + f_a(n * tau, (i + 1) * h)) / 2 * field[n][i]


def Wl(field, n, i, tau, h, f_a):
    return (f_a(n * tau, i * h) + f_a(n * tau, (i - 1) * h)) / 2 * field[n][i - 1]


def calculate_cell(field, n, i, tau, h, Wr, Wl, f_a):
    field[n][i] = field[n - 1][i] - tau / h * (Wr(field, n - 1, i, tau, h, f_a) - Wl(field, n - 1, i, tau, h, f_a))


class Solver:
    def __init__(self, T, tau, X, h, Wl, Wr, f_a, f_beg):
        self.tau, self.h = tau, h
        self.N_t, self.N_x = int(T / tau), int(X / h)
        self.Wl, self.Wr = Wl, Wr
        self.f_a = f_a
        self.Field = np.zeros((self.N_t, self.N_x), dtype=np.float64)
        self.Field[0] = [f_beg(h * i) for i in range(self.N_x)]

    def calculate(self):
        for n in range(1, self.N_t):
            for i in range(self.N_x):
                calculate_cell(self.Field, n, i, self.tau, self.h, Wr, Wl, self.f_a)

    def get_field(self):
        return self.Field.copy()

    def get_courant(self, t, x):
        return self.f_a(t, x) * self.tau / self.h


def get_acc_field(T, tau, X, h, f_a, f_beg):
    N_t, N_x = int(T / tau), int(X / h)
    Field = np.zeros((N_t, N_x), dtype=np.float64)
    Field[0] = [f_beg(h * i) for i in range(N_x)]
    for n in range(1, N_t):
        for i in range(N_x):
            Field[n][i] = f_beg((i*h - n*tau*f_a(n * tau, i * h)) % X)
    return Field