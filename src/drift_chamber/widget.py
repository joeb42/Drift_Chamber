import tkinter as tk
from tkinter import ttk

import matplotlib

matplotlib.use("TkAgg")
import copy

import matplotlib.animation as animation
import matplotlib.cm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure

from .driftchamber import DriftChamber
from .muon import MissedDetector, Muon


class Widget:
    def __init__(self):
        self.fig = Figure(figsize=(6, 4), dpi=150)
        self.sp = self.fig.add_subplot(111)
        self.sp.set_xlabel("y coordinate (cm)")
        self.sp.set_ylabel("z coordinate (cm)")
        self.sp.set_ylim(0, 30)
        self.sp.set_xlim(0, 50)
        self.root = tk.Tk()
        self.root.title("Drift Chamber")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().grid(column=0, row=3)
        spacings_label = tk.Label(self.root, text="Select a grid spacing")
        spacings_label.grid(row=0, column=1)
        self.spacings = ttk.Combobox(
            self.root, values=["1", "0.5", "0.1"], state="readonly"
        )
        self.spacings.grid(row=1, column=1)
        spacing_button = tk.Button(
            self.root, text="Press to confirm spacing", command=self.spacing_press
        )
        spacing_button.grid(row=2, column=1)
        muon_button = tk.Button(
            self.root, text="Press to generate random muon track", command=self.gen_muon
        )
        muon_button.grid(row=3, column=1)
        ani_button = tk.Button(
            self.root, text="Press to start animation", command=self.ani
        )
        ani_button.grid(row=4, column=1)
        stop_button = tk.Button(
            self.root, text="Press to stop animation", command=self.stop_press
        )
        stop_button.grid(row=5, column=1)
        clear_button = tk.Button(
            self.root, text="Press to clear axis", command=self.clear_press
        )
        clear_button.grid(row=6, column=1)
        tk.mainloop()

    def stop_press(self):
        self.stop = True

    def clear_press(self):
        self.chamber.grid *= 0.0
        self.fig.clf()
        self.sp = self.fig.add_subplot(111)
        self.sp.set_xlabel("y (cm)")
        self.sp.set_ylabel("z (cm)")
        self.sp.set_ylim(0, 30)
        self.sp.set_xlim(0, 50)
        self.canvas.draw()

    def spacing_press(self):
        self.chamber = DriftChamber(spacing=float(self.spacings.get()))

    def gen_muon(self):
        # Copy the jet cmap from matplotlib
        my_cmap = copy.copy(matplotlib.cm.get_cmap("jet"))
        # Set bad (zero value) pixels to dark blue
        my_cmap.set_bad((0, 0, 0.5))
        while True:
            try:
                m = Muon()
                self.chamber.initial_ion(m)
                self.im = self.sp.imshow(
                    self.chamber.grid,
                    cmap=my_cmap,
                    origin="lower",
                    extent=[0, 50, 0, 30],
                    norm=LogNorm(vmin=1e-2, vmax=100),
                )
                self.fig.colorbar(self.im, label=r"Charge density ($em^{-2}$)")
                self.canvas.draw()
                return
            except MissedDetector:
                continue

    def __animate(self, i):
        self.chamber.drift_diff()
        self.im.set_array(self.chamber.grid)
        self.fig.suptitle(f"t = {round(i * self.chamber.tau, 6)} seconds")
        return [self.im]

    def __gen(self):
        i = 0
        while not self.stop:
            i += 1
            yield i

    def ani(self):
        self.stop = False
        self.an = animation.FuncAnimation(
            self.fig, self.__animate, interval=15, frames=self.__gen, repeat=False
        )
        self.canvas.draw()
