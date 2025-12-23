from NM_interaction import InputApp
from NM_interaction import review_steel_strains, plot_M_contributions, plot_N_contributions, test_positive_rectangular_major_plot, test_positive_circular_major_plot

if __name__ == "__main__":
    import tkinter as tk
    root = tk.Tk()
    app = InputApp(root)
    root.mainloop()
    root.destroy()