import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Button

def mesh_parser(mesh_file):

    with open(mesh_file,'r') as file:
        line = file.readline()
        line = file.readline()
        line.strip()

        current_string = ""
        coord = []
        i = 0

        for c in line:
            i += 1
            if c == ',':
                coord.append(float(current_string))
                current_string = ""
            elif i == len(line):
                current_string += c
                coord.append(float(current_string))
            else:
                current_string += c

    return coord

    #print(str(coord[:]))

def output_parser(output_file):

    with open(output_file,'r') as file:
        line = file.readline()
        
        u = [[],[]]
        t = []
        i = 1

        for line in file:
            line.strip()

            u.append([])
            flag = True
            j = 0
            current_string = ""

            for c in line:
                j += 1
                if c == ',' and flag is False:
                    u[i].append(float(current_string))
                    current_string = ""
                elif c == ',' and flag is True:
                    t.append(float(current_string))
                    flag = False
                    current_string = ""
                elif j == len(line):
                    current_string += c
                    u[i].append(float(current_string))
                else:
                    current_string += c
            i += 1
    return t, u

def main():
    coord = mesh_parser("mesh.csv")
    t, u = output_parser("output.csv")

    x = np.array(coord)
    u = [np.array(ui) for ui in u]

    # --- Setup figura ---
    fig, ax = plt.subplots(figsize=(10, 5))
    line, = ax.plot([], [], lw=2, color='royalblue')
    time_text = ax.text(0.02, 0.92, '', transform=ax.transAxes, fontsize=12)

    ax.set_xlim(min(x), max(x))
    ax.set_ylim(min(v for frame in u for v in frame),
                max(v for frame in u for v in frame))
    ax.set_xlabel("x")
    ax.set_ylabel("u")
    ax.set_title("Soluzione 1D")
    ax.grid(True)

    # --- Init: pulisce il grafico prima di partire ---
    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    # --- Update: chiamata ad ogni frame con l'indice del frame corrente ---
    def update(frame):
        line.set_data(x, u[frame])
        time_text.set_text(f't = {t[frame]:.4f}')
        return line, time_text

    # --- Animazione ---
    ani = animation.FuncAnimation(
        fig,
        update,
        frames=len(t),
        init_func=init,
        interval=50,   # ms tra un frame e l'altro — abbassa per velocizzare
        blit=True
    )

    # --- Bottone Pause/Play ---
    paused = False  # ← variabile locale in main()

    plt.subplots_adjust(bottom=0.15)
    ax_button = plt.axes([0.45, 0.02, 0.12, 0.05])
    btn = Button(ax_button, 'Pause')

    def toggle(event):
        nonlocal paused  # ← accede alla variabile di main(), non globale
        if paused:
            ani.resume()
            btn.label.set_text('Pause')
        else:
            ani.pause()
            btn.label.set_text('Play')
        paused = not paused
        plt.draw()

    btn.on_clicked(toggle)

    plt.tight_layout()
    plt.show()




if __name__ == "__main__":
    print("""
    =========================================
      Visualizzatore Risultati Solver 1D
    =========================================
    """)
    main()