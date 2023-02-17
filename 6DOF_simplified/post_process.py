import matplotlib.pyplot as plt

def read_data(file_path):
    t = []
    x = []
    y = []
    z = []
    rx = []
    ry = []
    rz = []
    vx = []
    vy = []
    vz = []
    u = []
    with open(file_path, "r") as file:
        for line in file:
            data = line.strip().split()
            t.append(float(data[0]))
            x.append(float(data[1]))
            y.append(float(data[2]))
            z.append(float(data[3]))
            rx.append(float(data[4]))
            ry.append(float(data[5]))
            rz.append(float(data[6]))
            vx.append(float(data[7]))
            vy.append(float(data[8]))
            vz.append(float(data[9]))
            u.append(float(data[10]))
    return t, x, y, z, rx, ry, rz, vx, vy, vz, u;

file_path = "cmake-build-release/simulation_results.txt"
t, x, y, z, rx, ry, rz, vx, vy, vz, u = read_data(file_path)


fig, axs = plt.subplots(2, 2)

axs[0, 0].plot(t, vx, label="Vx")
axs[0, 0].plot(t, vy, label="Vy")
axs[0, 0].plot(t, vz, label="Vz")
axs[0, 0].set_title("Velocities")
axs[0, 0].legend(loc="best")

axs[0, 1].plot(t, rx, label="rx")
axs[0, 1].plot(t, ry, label="ry")
axs[0, 1].plot(t, rz, label="rz")
axs[0, 1].set_title("Angles")
axs[0, 1].legend(loc="best")

axs[1, 0].plot(x, z, label="z")
axs[1, 0].set_title("altitude")
axs[1, 0].legend(loc="best")

axs[1, 1].plot(x, y, label="position")
axs[1, 1].set_title("position")
axs[1, 1].legend(loc="best")


plt.tight_layout()
plt.show()

plt.plot(t, u, label="input")
plt.legend()
plt.show()