import os
import tkinter
import customtkinter

from PIL import Image

customtkinter.set_appearance_mode("System")         # Modes: "System" (standard), "Dark", "Light"
customtkinter.set_default_color_theme("dark-blue")  # Themes: "blue" (standard), "green", "dark-blue"

class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        # configure window
        width = 800
        height = 200

        self.title("SeisFAT3D - Seismic First Arrival Toolkit 3D")
        self.geometry(f"{width}x{height}+200+100")  
        self.resizable(False, False)

        compile_buttom = customtkinter.CTkButton(self, height = 50, width = 180, text = "Compile", font=customtkinter.CTkFont(size=20, weight="normal"))
        compile_buttom.place(x = 10, y = 10)

        modeling_buttom = customtkinter.CTkButton(self, height = 50, width = 180, text = "Modeling", font=customtkinter.CTkFont(size=20, weight="normal"))
        modeling_buttom.place(x = 210, y = 10)

        inverion_buttom = customtkinter.CTkButton(self, height = 50, width = 180, text = "Inversion", font=customtkinter.CTkFont(size=20, weight="normal"))
        inverion_buttom.place(x = 410, y = 10)

        migration_buttom = customtkinter.CTkButton(self, height = 50, width = 180, text = "Migration", font=customtkinter.CTkFont(size=20, weight="normal"))
        migration_buttom.place(x = 610, y = 10)


        tools_buttom = customtkinter.CTkButton(self, height = 30, width = 180, text = "Tools", font=customtkinter.CTkFont(size=20, weight="normal"))
        tools_buttom.place(x = 10, y = 160)

        help_buttom = customtkinter.CTkButton(self, height = 30, width = 180, text = "Help", font=customtkinter.CTkFont(size=20, weight="normal"))
        help_buttom.place(x = 210, y = 160)


        appearance_mode_label = customtkinter.CTkLabel(self, text="Appearance Mode", font = customtkinter.CTkFont(size=18))
        appearance_mode_label.place(x = 630, y = 120)
        appearance_mode_optionemenu = customtkinter.CTkOptionMenu(self, height = 30, width = 180, font = customtkinter.CTkFont(size=18), values = ["System", "Light", "Dark"], command = customtkinter.set_appearance_mode)
        appearance_mode_optionemenu.place(x = 610, y = 160)
        


        # image_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "imgs/")
        # self.logo_image = customtkinter.CTkImage(Image.open(os.path.join(image_path, "logo.png")), size=(180, 130))

        # self.logo_frame = customtkinter.CTkFrame(self.menu_frame, height = 130, width = 180, corner_radius=0) 
        # self.logo_frame.place(x = 20, y = 20)

        # self.img_frame_label = customtkinter.CTkLabel(self.logo_frame, text = " ", image=self.logo_image)
        # self.img_frame_label.place(x = 0, y = 0)

        # self.name_frame = customtkinter.CTkFrame(self.logo_frame, fg_color = ("#EEEEEE","#191A19"), height = 30, width = 180, corner_radius=0)    
        # self.name_frame.place(x = 0, y = 100)

        # self.name_label = customtkinter.CTkLabel(self.name_frame, text="HeadTomoGrid 3D", font=customtkinter.CTkFont(size=15))
        # self.name_label.place(relx = 0.5, rely = 0.5, anchor = "center")

#-----------------------------------------------------------------------------------------
# Application buttons
#-----------------------------------------------------------------------------------------

    #     self.model_buttom = customtkinter.CTkButton(self.menu_frame, height = 40, width = 180, text = "Model", font=customtkinter.CTkFont(size=20, weight="normal"), command=self.model_button_event)
    #     self.model_buttom.place(x = 20, rely = 200/self.height)

    #     self.geometry_buttom = customtkinter.CTkButton(self.menu_frame, height = 40, width = 180, text = "Geometry", font=customtkinter.CTkFont(size=20, weight="normal"))#, command=self.geometry_button_event)
    #     self.geometry_buttom.place(x = 20, rely = 260/self.height)

    #     self.eikonal_buttom = customtkinter.CTkButton(self.menu_frame, height = 40, width = 180, text = "Eikonal", font=customtkinter.CTkFont(size=20, weight="normal"))#, command=self.source_button_event)
    #     self.eikonal_buttom.place(x = 20, rely = 320/self.height)

    #     self.tomography_buttom = customtkinter.CTkButton(self.menu_frame, height = 40, width = 180, text = "Tomography", font=customtkinter.CTkFont(size=20, weight="normal"))#, command=self.tomography_button_event)
    #     self.tomography_buttom.place(x = 20, rely = 380/self.height)

        # self.appearance_mode_label = customtkinter.CTkLabel(self.menu_frame, text="Appearance Mode", font=customtkinter.CTkFont(size=17))
        # self.appearance_mode_label.place(relx = 0.5, y = 0, anchor = "center")
        # self.appearance_mode_optionemenu = customtkinter.CTkOptionMenu(self.menu_frame, height = 30, width = 180, font=customtkinter.CTkFont(size=15), values=["Light", "Dark", "System"],command=self.change_appearance_mode_event)
        # self.appearance_mode_optionemenu.place(relx = 0.5, y = 8, anchor = "center")
        # self.appearance_mode_optionemenu.set("System")

    #     self.main_frame = customtkinter.CTkFrame(self, width = 1240, height = 800, fg_color = "transparent")    

    # def change_appearance_mode_event(self, new_appearance_mode: str):
    #     customtkinter.set_appearance_mode(new_appearance_mode)

    # def main_frame_clear(self):
    #     for widget in self.main_frame.winfo_children():
    #         widget.destroy()

    #     self.main_frame.place_forget()
    #     self.main_frame.place(x = 240, y = 20)     

#-----------------------------------------------------------------------------------------
# Model -------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

    # def model_button_event(self):
    #     print("model_button click")
    #     self.main_frame_clear()

    #     self.model_dimension_frame = customtkinter.CTkFrame(self.main_frame, width = 400, height = 250, corner_radius = 5) 
    #     self.model_dimension_frame.place(x = 0, y = 0)    

    #     self.dimension_title_frame = customtkinter.CTkFrame(self.model_dimension_frame, width = 400, height = 50, corner_radius = 5)
    #     self.dimension_title_frame.place(x = 0, y = 0)    

    #     self.dimension_title_label = customtkinter.CTkLabel(self.dimension_title_frame, text="Dimensions", font=customtkinter.CTkFont(size=20))
    #     self.dimension_title_label.place(relx = 0.5, rely = 0.5, anchor= "center")

    #     # Model samples 

    #     self.model_samples_frame = customtkinter.CTkFrame(self.model_dimension_frame, width = 180, height = 180, corner_radius = 0)
    #     self.model_samples_frame.place(x = 10, y = 60) 

    #     self.model_samples_title_frame = customtkinter.CTkFrame(self.model_samples_frame, fg_color = ("#EEEEEE","#191A19"), width = 180, height = 30, corner_radius = 0)
    #     self.model_samples_title_frame.place(x = 0, y = 0) 

    #     self.model_samples_title = customtkinter.CTkLabel(self.model_samples_title_frame, text="Samples", font=customtkinter.CTkFont(size=15))
    #     self.model_samples_title.place(relx = 0.5, rely = 0.5, anchor = "center")

    #     self.Nx_label = customtkinter.CTkLabel(self.model_samples_frame, text="Nx:", font=customtkinter.CTkFont(size=15))
    #     self.Nx_label.place(x = 20, y = 50)
    #     self.Nx_entry = customtkinter.CTkEntry(self.model_samples_frame, width = 110, placeholder_text="x unit samples", font=customtkinter.CTkFont(size=15))
    #     self.Nx_entry.place(x = 50, y = 50)

    #     self.Ny_label = customtkinter.CTkLabel(self.model_samples_frame, text="Ny:", font=customtkinter.CTkFont(size=15))
    #     self.Ny_label.place(x = 20, y = 90)
    #     self.Ny_entry = customtkinter.CTkEntry(self.model_samples_frame, width = 110, placeholder_text="y unit samples", font=customtkinter.CTkFont(size=15))
    #     self.Ny_entry.place(x = 50, y = 90)

    #     self.Nz_label = customtkinter.CTkLabel(self.model_samples_frame, text="Nz:", font=customtkinter.CTkFont(size=15))
    #     self.Nz_label.place(x = 20, y = 130)
    #     self.Nz_entry = customtkinter.CTkEntry(self.model_samples_frame, width = 110, placeholder_text="z unit samples", font=customtkinter.CTkFont(size=15))
    #     self.Nz_entry.place(x = 50, y = 130)

    #     # Model spacing

    #     self.model_spacing_frame = customtkinter.CTkFrame(self.model_dimension_frame, width = 180, height = 180, corner_radius = 0)
    #     self.model_spacing_frame.place(x = 210, y = 60) 

    #     self.model_spacing_title_frame = customtkinter.CTkFrame(self.model_spacing_frame, fg_color = ("#EEEEEE","#191A19"), width = 180, height = 30, corner_radius = 0)
    #     self.model_spacing_title_frame.place(x = 0, y = 0) 

    #     self.model_spacing_title = customtkinter.CTkLabel(self.model_spacing_title_frame, text="Spacing", font=customtkinter.CTkFont(size=15))
    #     self.model_spacing_title.place(relx = 0.5, rely = 0.5, anchor = "center")

    #     self.dx_label = customtkinter.CTkLabel(self.model_spacing_frame, text="dx:", font=customtkinter.CTkFont(size=15))
    #     self.dx_label.place(x = 20, y = 50)
    #     self.dx_entry = customtkinter.CTkEntry(self.model_spacing_frame, width = 110, placeholder_text="x spacing [m]", font=customtkinter.CTkFont(size=15))
    #     self.dx_entry.place(x = 50, y = 50)

    #     self.dy_label = customtkinter.CTkLabel(self.model_spacing_frame, text="dy:", font=customtkinter.CTkFont(size=15))
    #     self.dy_label.place(x = 20, y = 90)
    #     self.dy_entry = customtkinter.CTkEntry(self.model_spacing_frame, width = 110, placeholder_text="y spacing [m]", font=customtkinter.CTkFont(size=15))
    #     self.dy_entry.place(x = 50, y = 90)

    #     self.dz_label = customtkinter.CTkLabel(self.model_spacing_frame, text="Nz:", font=customtkinter.CTkFont(size=15))
    #     self.dz_label.place(x = 20, y = 130)
    #     self.dz_entry = customtkinter.CTkEntry(self.model_spacing_frame, width = 110, placeholder_text="z spacing [m]", font=customtkinter.CTkFont(size=15))
    #     self.dz_entry.place(x = 50, y = 130)

    #     # Model visualization

    #     self.model_visualization_frame = customtkinter.CTkFrame(self.main_frame, width = 400, height = 570, corner_radius = 5) 
    #     self.model_visualization_frame.place(x = 0, y = 270)    

    #     self.model_visualization_title_frame = customtkinter.CTkFrame(self.model_visualization_frame, width = 400, height = 50, corner_radius = 5)
    #     self.model_visualization_title_frame.place(x = 0, y = 0)    

    #     self.model_visualization_title_label = customtkinter.CTkLabel(self.model_visualization_title_frame, text="Visualization", font=customtkinter.CTkFont(size=20))
    #     self.model_visualization_title_label.place(relx = 0.5, rely = 0.5, anchor= "center")

    #     # slices

    #     self.model_visualization_slices_frame = customtkinter.CTkFrame(self.model_visualization_frame, width = 180, height = 180, corner_radius = 0)
    #     self.model_visualization_slices_frame.place(x = 10, y = 60) 

    #     self.model_visualization_slices_title_frame = customtkinter.CTkFrame(self.model_visualization_slices_frame, fg_color = ("#EEEEEE","#191A19"), width = 180, height = 30, corner_radius = 0)
    #     self.model_visualization_slices_title_frame.place(x = 0, y = 0) 

    #     self.model_visualization_slices_title = customtkinter.CTkLabel(self.model_visualization_slices_title_frame, text="Slices", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_slices_title.place(relx = 0.5, rely = 0.5, anchor = "center")

    #     self.model_sliceZY_label = customtkinter.CTkLabel(self.model_visualization_slices_frame, text="ZY:", font=customtkinter.CTkFont(size=15))
    #     self.model_sliceZY_label.place(x = 20, y = 50)
    #     self.model_sliceZY_entry = customtkinter.CTkEntry(self.model_visualization_slices_frame, width = 110, placeholder_text="x position [m]", font=customtkinter.CTkFont(size=15))
    #     self.model_sliceZY_entry.place(x = 50, y = 50)

    #     self.model_sliceZX_label = customtkinter.CTkLabel(self.model_visualization_slices_frame, text="ZX:", font=customtkinter.CTkFont(size=15))
    #     self.model_sliceZX_label.place(x = 20, y = 90)
    #     self.model_sliceZX_entry = customtkinter.CTkEntry(self.model_visualization_slices_frame, width = 110, placeholder_text="y position [m]", font=customtkinter.CTkFont(size=15))
    #     self.model_sliceZX_entry.place(x = 50, y = 90)

    #     self.model_sliceXY_label = customtkinter.CTkLabel(self.model_visualization_slices_frame, text="XY:", font=customtkinter.CTkFont(size=15))
    #     self.model_sliceXY_label.place(x = 20, y = 130)
    #     self.model_sliceXY_entry = customtkinter.CTkEntry(self.model_visualization_slices_frame, width = 110, placeholder_text="z position [m]", font=customtkinter.CTkFont(size=15))
    #     self.model_sliceXY_entry.place(x = 50, y = 130)

    #     # geometry

    #     self.model_visualization_geometry_frame = customtkinter.CTkFrame(self.model_visualization_frame, width = 180, height = 180, corner_radius = 0)
    #     self.model_visualization_geometry_frame.place(x = 210, y = 60) 

    #     self.model_visualization_geometry_title_frame = customtkinter.CTkFrame(self.model_visualization_geometry_frame, fg_color = ("#EEEEEE","#191A19"), width = 180, height = 30, corner_radius = 0)
    #     self.model_visualization_geometry_title_frame.place(x = 0, y = 0) 

    #     self.model_visualization_geometry_title = customtkinter.CTkLabel(self.model_visualization_geometry_title_frame, text="Geometry", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_geometry_title.place(relx = 0.5, rely = 0.5, anchor = "center")

    #     self.model_visualization_shots_check = customtkinter.CTkCheckBox(self.model_visualization_geometry_frame, text=" View shots", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_shots_check.place(relx = 0.5, y = 65, anchor = "center")

    #     self.model_visualization_nodes_check = customtkinter.CTkCheckBox(self.model_visualization_geometry_frame, text=" View nodes", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_nodes_check.place(relx = 0.5, y = 105, anchor = "center")

    #     self.model_visualization_cmps_check = customtkinter.CTkCheckBox(self.model_visualization_geometry_frame, text=" View CMPs", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_cmps_check.configure(state = "disabled")
    #     self.model_visualization_cmps_check.place(relx = 0.5, y = 145, anchor = "center")

    #     # Filenames

    #     self.model_visualization_filenames_frame = customtkinter.CTkFrame(self.model_visualization_frame, width = 380, height = 220, corner_radius = 0)
    #     self.model_visualization_filenames_frame.place(x = 10, y = 250)

    #     self.model_visualization_filenames_title_frame = customtkinter.CTkFrame(self.model_visualization_filenames_frame, fg_color = ("#EEEEEE","#191A19"), width = 380, height = 30, corner_radius = 0)
    #     self.model_visualization_filenames_title_frame.place(x = 0, y = 0)

    #     self.model_visualization_filenames_title_label = customtkinter.CTkLabel(self.model_visualization_filenames_title_frame, text="Filenames", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_filenames_title_label.place(relx = 0.5, rely = 0.5, anchor = "center")

    #     self.model_visualization_filenames_model_label = customtkinter.CTkLabel(self.model_visualization_filenames_frame, text="Model:", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_filenames_model_label.place(x = 20, y = 50)
    #     self.model_visualization_filenames_model_entry = customtkinter.CTkEntry(self.model_visualization_filenames_frame, width = 285, placeholder_text="Model location [.bin] (ord = 'F')", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_filenames_model_entry.place(x = 75, y = 50)

    #     self.model_visualization_filenames_shots_label = customtkinter.CTkLabel(self.model_visualization_filenames_frame, text="Shots:", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_filenames_shots_label.place(x = 20, y = 90)
    #     self.model_visualization_filenames_shots_entry = customtkinter.CTkEntry(self.model_visualization_filenames_frame, width = 285, placeholder_text="Shots location [xyz.txt]", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_filenames_shots_entry.place(x = 75, y = 90)
    #     self.model_visualization_filenames_shots_entry.configure(state = "disabled")

    #     self.model_visualization_filenames_nodes_label = customtkinter.CTkLabel(self.model_visualization_filenames_frame, text="Nodes:", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_filenames_nodes_label.place(x = 20, y = 130)
    #     self.model_visualization_filenames_nodes_entry = customtkinter.CTkEntry(self.model_visualization_filenames_frame, width = 285, placeholder_text="Nodes location [xyz.txt]", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_filenames_nodes_entry.place(x = 75, y = 130)
    #     self.model_visualization_filenames_nodes_entry.configure(state = "disabled")

    #     self.model_visualization_filenames_figure_label = customtkinter.CTkLabel(self.model_visualization_filenames_frame, text="Figure:", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_filenames_figure_label.place(x = 20, y = 170)
    #     self.model_visualization_filenames_figure_entry = customtkinter.CTkEntry(self.model_visualization_filenames_frame, width = 285, placeholder_text="Figure location [.png]", font=customtkinter.CTkFont(size=15))
    #     self.model_visualization_filenames_figure_entry.place(x = 75, y = 170)

    #     # Apply button

    #     self.model_visualization_apply_buttom = customtkinter.CTkButton(self.model_visualization_frame, height = 40, width = 120, text = "Show", font=customtkinter.CTkFont(size=20, weight="normal"))
    #     self.model_visualization_apply_buttom.place(x = 270, y = 480)

    #     # Building a simple model

    #     self.model_building_frame = customtkinter.CTkFrame(self.main_frame, width = 820, height = 390, corner_radius = 5) 
    #     self.model_building_frame.place(x = 420, y = 0)    

    #     self.building_title_frame = customtkinter.CTkFrame(self.model_building_frame, width = 820, height = 50, corner_radius = 5)
    #     self.building_title_frame.place(x = 0, y = 0)    

    #     self.building_title_label = customtkinter.CTkLabel(self.building_title_frame, text="Simple Model Building", font=customtkinter.CTkFont(size=20))
    #     self.building_title_label.place(relx = 0.5, rely = 0.5, anchor = "center")

    #     # layer cake

    #     self.model_building_layercake_frame = customtkinter.CTkFrame(self.model_building_frame, width = 300, height = 320, corner_radius = 0) 
    #     self.model_building_layercake_frame.place(x = 10, y = 60)    

    #     self.model_building_layercake_title_frame = customtkinter.CTkFrame(self.model_building_layercake_frame, width = 300, height = 30, fg_color = ("#EEEEEE","#191A19"), corner_radius = 0)
    #     self.model_building_layercake_title_frame.place(x = 0, y = 0)    

    #     self.model_building_layercake_title_label = customtkinter.CTkLabel(self.model_building_layercake_title_frame, text="Layer cake", font=customtkinter.CTkFont(size=15))
    #     self.model_building_layercake_title_label.place(relx = 0.5, rely = 0.5, anchor= "center")






    #     self.model_building_layercake_velocity_label = customtkinter.CTkLabel(self.model_building_layercake_frame, text="Vp:", font=customtkinter.CTkFont(size=15)) 
    #     self.model_building_layercake_velocity_label.place(x = 45, y = 200)
    #     self.model_building_layercake_velocity_entry = customtkinter.CTkEntry(self.model_building_layercake_frame, width = 200, placeholder_text="Velocity array [m/s]", font=customtkinter.CTkFont(size=15)) 
    #     self.model_building_layercake_velocity_entry.place(x = 80, y = 200)

    #     self.model_building_layercake_depth_label = customtkinter.CTkLabel(self.model_building_layercake_frame, text="Depth:", font=customtkinter.CTkFont(size=15)) 
    #     self.model_building_layercake_depth_label.place(x = 25, y = 240)
    #     self.model_building_layercake_depth_entry = customtkinter.CTkEntry(self.model_building_layercake_frame, width = 200, placeholder_text="Interface array [m]", font=customtkinter.CTkFont(size=15)) 
    #     self.model_building_layercake_depth_entry.place(x = 80, y = 240)

    #     self.model_building_layercake_output_label = customtkinter.CTkLabel(self.model_building_layercake_frame, text="Output:", font=customtkinter.CTkFont(size=15)) 
    #     self.model_building_layercake_output_label.place(x = 20, y = 280)
    #     self.model_building_layercake_output_entry = customtkinter.CTkEntry(self.model_building_layercake_frame, width = 200, placeholder_text="Save location [.bin] (ord = 'F')", font=customtkinter.CTkFont(size=15)) 
    #     self.model_building_layercake_output_entry.place(x = 80, y = 280)

    #     # linear gradient

    #     self.model_building_gradient_frame = customtkinter.CTkFrame(self.model_building_frame, width = 300, height = 320, corner_radius = 0) 
    #     self.model_building_gradient_frame.place(x = 510, y = 60)    

    #     self.model_building_gradient_title_frame = customtkinter.CTkFrame(self.model_building_gradient_frame, width = 300, height = 30, fg_color = ("#EEEEEE","#191A19"), corner_radius = 0)
    #     self.model_building_gradient_title_frame.place(x = 0, y = 0)    

    #     self.model_building_gradient_title_label = customtkinter.CTkLabel(self.model_building_gradient_title_frame, text="Linear gradient", font=customtkinter.CTkFont(size=15))
    #     self.model_building_gradient_title_label.place(relx = 0.5, rely = 0.5, anchor= "center")


    #     self.model_building_gradient_input_label = customtkinter.CTkLabel(self.model_building_gradient_frame, text="Input:", font=customtkinter.CTkFont(size=15)) 
    #     self.model_building_gradient_input_label.place(x = 30, y = 240)
    #     self.model_building_gradient_input_entry = customtkinter.CTkEntry(self.model_building_gradient_frame, width = 200, placeholder_text="Location [.bin] (ord = 'F')", font=customtkinter.CTkFont(size=15)) 
    #     self.model_building_gradient_input_entry.place(x = 80, y = 240)

    #     self.model_building_gradient_output_label = customtkinter.CTkLabel(self.model_building_gradient_frame, text="Output:", font=customtkinter.CTkFont(size=15)) 
    #     self.model_building_gradient_output_label.place(x = 20, y = 280)
    #     self.model_building_gradient_output_entry = customtkinter.CTkEntry(self.model_building_gradient_frame, width = 200, placeholder_text="Location [.bin] (ord = 'F')", font=customtkinter.CTkFont(size=15)) 
    #     self.model_building_gradient_output_entry.place(x = 80, y = 280)

    #     # model type 

    #     self.model_building_type_frame = customtkinter.CTkFrame(self.model_building_frame, width = 180, height = 145, corner_radius = 0) 
    #     self.model_building_type_frame.place(x = 320, y = 60)    

    #     self.model_building_type_title_frame = customtkinter.CTkFrame(self.model_building_type_frame, width = 180, height = 30, fg_color = ("#EEEEEE","#191A19"), corner_radius = 0)
    #     self.model_building_type_title_frame.place(x = 0, y = 0)    

    #     self.model_building_type_title_label = customtkinter.CTkLabel(self.model_building_type_title_frame, text="Model type", font=customtkinter.CTkFont(size=15))
    #     self.model_building_type_title_label.place(relx = 0.5, rely = 0.5, anchor= "center")

    #     self.model_building_layercake_check = customtkinter.CTkCheckBox(self.model_building_type_frame, text=" Layer cake", font=customtkinter.CTkFont(size=15))
    #     self.model_building_layercake_check.place(x = 20, y = 50)

    #     self.model_building_gradient_check = customtkinter.CTkCheckBox(self.model_building_type_frame, text=" Linear gradient", font=customtkinter.CTkFont(size=15))
    #     self.model_building_gradient_check.place(x = 20, y = 100)

    #     # build buttom

    #     self.model_building_apply_buttom = customtkinter.CTkButton(self.model_building_frame, height = 40, width = 180, text = "Build", font=customtkinter.CTkFont(size=20, weight="normal"))
    #     self.model_building_apply_buttom.place(x = 320, y = 215)

    #     # Model transformation

    #     self.model_transformation_frame = customtkinter.CTkFrame(self.main_frame, width = 820, height = 390, corner_radius = 5) 
    #     self.model_transformation_frame.place(x = 420, y = 410)    

    #     self.transformation_title_frame = customtkinter.CTkFrame(self.model_transformation_frame, width = 820, height = 50, corner_radius = 5)
    #     self.transformation_title_frame.place(x = 0, y = 0)    

    #     self.transformation_title_label = customtkinter.CTkLabel(self.transformation_title_frame, text="Model Transformation", font=customtkinter.CTkFont(size=20))
    #     self.transformation_title_label.place(relx = 0.5, rely = 0.5, anchor= "center")








# #-----------------------------------------------------------------------------------------
# # Geometry -------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def geometry_button_event(self):
#         print("geometry_button click")
#         self.main_frame_clear()


# #-----------------------------------------------------------------------------------------
# # Source ---------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def source_button_event(self):
#         print("source_button click")
#         self.main_frame_clear()

# #-----------------------------------------------------------------------------------------
# # Acoustic -------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def acoustic_button_event(self):
#         print("acoustic_button click")
#         self.main_frame_clear()

# #-----------------------------------------------------------------------------------------
# # Auto-picking ---------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def autoPick_button_event(self):
#         print("autoPick_button click")
#         self.main_frame_clear()

# #-----------------------------------------------------------------------------------------
# # Eikonal --------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def eikonal_button_event(self):
#         print("eikonal_button click")
#         self.main_frame_clear()

# #-----------------------------------------------------------------------------------------
# # Tomography -----------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------

#     def tomography_button_event(self):
#         print("tomography_button click")
#         self.main_frame_clear()


#-----------------------------------------------------------------------------------------
# Main function --------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

if __name__ == "__main__":
    app = App()
    app.mainloop()
