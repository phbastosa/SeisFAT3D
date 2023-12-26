import customtkinter as ctk

class Level():
    def __init__(self, window, text, ypos):
        
        self.button = []
        
        self.frame = ctk.CTkFrame(window, height = 80, width = 800)
        self.frame.place(x = 100, y = ypos)

        self.title_frame = ctk.CTkFrame(self.frame, height = 30, width = 300)
        self.title_frame.place(x = 0, y = 0)

        self.title_label = ctk.CTkLabel(self.title_frame, text = text, font = ("Helvetica", 20))
        self.title_label.place(x = 10, y = 2)

class Application():
    def __init__(self):
    
        self.set_main_window()
        self.set_main_frames()
        self.set_main_buttons()
        self.set_initial_levels()

        self.root.mainloop()

    def set_main_window(self):

        self.root = ctk.CTk()
        self.root.title("SeisFAT3D - Seismic First Arrival Toolkit 3D")

        self.myFont = ("Helvetica", 20)

        self.window_width = 1000
        self.window_height = self.root.winfo_screenheight()

        self.root.minsize(self.window_width, self.window_height)
        self.root.maxsize(self.window_width, self.window_height)

        center_x = int(0.5*self.root.winfo_screenwidth() - 0.5*self.window_width)

        self.root.geometry(f'{self.window_width}x{self.window_height}+{center_x}+0')

    def set_main_frames(self):

        self.upper_frame = ctk.CTkFrame(self.root, height = 40, width = self.window_width)
        self.upper_frame.place(x = 0, y = 0)

        self.lower_frame = ctk.CTkFrame(self.root, height = 100, width = self.window_width)
        self.lower_frame.place(x = 0, y = self.window_height - 105)

    def set_main_buttons(self):
        
        self.reset_button = ctk.CTkButton(self.upper_frame, height = 30, width = 120, text = "Reset", font = self.myFont)
        self.reset_button.place(x = 5, y = 5)
        
        self.delete_button = ctk.CTkButton(self.upper_frame, height = 30, width = 120, text = "Delete", font = self.myFont)
        self.delete_button.place(x = 135, y = 5)

        self.add_button = ctk.CTkButton(self.upper_frame, height = 30, width = 120, text = "Add", font = self.myFont)
        self.add_button.place(x = 265, y = 5)

        self.load_button = ctk.CTkButton(self.lower_frame, height = 30, width = 120, text = "Load", font = self.myFont)
        self.load_button.place(x = 5, y = 5)
        
        self.save_button = ctk.CTkButton(self.lower_frame, height = 30, width = 120, text = "Save", font = self.myFont)
        self.save_button.place(x = 135, y = 5)

        self.run_button = ctk.CTkButton(self.lower_frame, height = 30, width = 120, text = "Run", font = self.myFont)
        self.run_button.place(x = self.window_width - 125, y = 5)

    def set_initial_levels(self):

        self.level_positions = [100, 200, 300, 400, 500, 600]
        self.level_names = ["Compilation", "Model & Geometry", "Synthetic Seismogram",
                            "First Arrival Modeling", "Tomography Inversion", "Kirchhoff Migration"]
        
        self.levels = []
        for i in range(len(self.level_positions)):
            self.levels.append(Level(self.root, f"Level {i+1}: {self.level_names[i]}", self.level_positions[i])) 

        self.compile_button = ctk.CTkButton(self.levels[0].frame, height = 30, width = 120, text = "Compile", font = self.myFont)
        self.compile_button.place(relx = 0.5, y = 55, anchor = "center")

        self.model_button = ctk.CTkButton(self.levels[1].frame, height = 30, width = 120, text = "Model", font = self.myFont)
        self.model_button.place(x = 400 - 65, y = 55, anchor = "center") 

        self.geometry_button = ctk.CTkButton(self.levels[1].frame, height = 30, width = 120, text = "Geometry", font = self.myFont)
        self.geometry_button.place(x = 400 + 65, y = 55, anchor = "center") 

    # def set_multiple_buttons(self, index):

    #     self.levels[index].frame.update_idletasks()

    #     ws = 50
    #     fw = 800
    #     bw = 120
    #     tw = sum(bw + ws for _ in range(len(self.levels[index].button)))

    #     init_relx = 0.5*(1.0 - (tw - 0.3*ws) / fw)

    #     for i in range(len(self.levels[1].button)):

    #         rx = init_relx + i*(bw + ws) / fw 

    #         self.levels[1].button[i].place(relx = rx, rely = 0.4, anchor = "w")

    def add_process(self):
        # new window     
        pass        

    def delete_process(self):
        # new window
        pass    

    def reset_process(self):
        # just a screen refresh
        pass    






if __name__ == "__main__":
    Application()