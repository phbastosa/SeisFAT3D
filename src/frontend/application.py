import tkinter as tk
import tkinter.font as font

root = tk.Tk()

root.title("SeisFAT3D - Seismic First Arrival Toolkit 3D")

window_width = 940
window_height = 500

myFont = font.Font(family='Helvetica', size = 15)

root.minsize(window_width, window_height)

# get the screen dimension
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

# find the center point
center_x = int(0.25*screen_width - 0.5*window_width)
center_y = int(0.25*screen_height - 0.5*window_height)

# set the position of the window to the center of the screen
root.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')

upper_frame = tk.Frame(root, height = 45, width = window_width, highlightthickness = True)
upper_frame.place(x = 0, y = 0, relwidth = 1.0)

load_button = tk.Button(upper_frame, height = 1, width = 10, bg = "#A9A9A9", text = "Load", justify = "center")
load_button.pack(side = "left", padx = 2, pady = 2)
load_button['font'] = myFont

save_button = tk.Button(upper_frame, height = 1, width = 10, bg = "#A9A9A9", text = "Save", justify = "center")
save_button.pack(side = "left", padx = 2, pady = 2)
save_button['font'] = myFont

reset_button = tk.Button(upper_frame, height = 1, width = 10, bg = "#A9A9A9", text = "Reset", justify = "center")
reset_button.pack(side = "right", padx = 2, pady = 2)
reset_button['font'] = myFont

del_button = tk.Button(upper_frame, height = 1, width = 10, bg = "#A9A9A9", text = "Delete", justify = "center")
del_button.pack(side = "right", padx = 2, pady = 2)
del_button['font'] = myFont

add_button = tk.Button(upper_frame, height = 1, width = 10, bg = "#A9A9A9", text = "Add", justify = "center")
add_button.pack(side = "right", padx = 2, pady = 2)
add_button['font'] = myFont


level0_frame = tk.LabelFrame(root, height = 60, text = "Level 0: Compilation")
level0_frame.place(relx = 0.2, y = 100, relwidth = 0.6)

compile_button = tk.Button(level0_frame, height = 1, width = 10, bg = "#A9A9A9", text = "Compile", justify = "center")
compile_button.place(anchor = "center", relx = 0.5, rely = 0.35)
compile_button['font'] = myFont

level1_frame = tk.LabelFrame(root, height = 30, text = "Level 1: Model & Geometry")
level1_frame.place(relx = 0.2, y = 180, relwidth = 0.6)

level2_frame = tk.LabelFrame(root, height = 30, text = "Level 2: Wavefield modeling")
level2_frame.place(relx = 0.2, y = 230, relwidth = 0.6)

level3_frame = tk.LabelFrame(root, height = 30, text = "Level 3: First Arrivals modeling")
level3_frame.place(relx = 0.2, y = 280, relwidth = 0.6)

level4_frame = tk.LabelFrame(root, height = 30, text = "Level 4: First Arrival tomography")
level4_frame.place(relx = 0.2, y = 330, relwidth = 0.6)

level5_frame = tk.LabelFrame(root, height = 30, text = "Level 5: Kirchhoff migration")
level5_frame.place(relx = 0.2, y = 380, relwidth = 0.6)

root.mainloop()