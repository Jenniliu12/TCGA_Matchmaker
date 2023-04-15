import tkinter as tk
from tkinter import ttk
from tkinter import *
from tkinter import filedialog
import pandas as pd
import numpy as np
from functools import partial
import match_computation as m

# START WINDOW EVENT LOOP:
window = tk.Tk()



input_ref_profile_var = tk.StringVar()
input_tcga_profile_var = tk.StringVar()
missing_values_var = tk.StringVar()
distance_metrics_var = tk.StringVar()
show_output_var = tk.StringVar()

global missing_values
def enter_entry():
    missing_values = missing_values_var.get()
    distance_metrics = distance_metrics_var.get()
    show_output = show_output_var.get()

    # RESET
    #missing_values_var.set("") 
    #distance_metrics_var.set("")
    #show_output_var.set("")


def enter_profiles():
    input_ref_profile = input_ref_profile_var.get()
    input_tcga_profile = input_tcga_profile_var.get()

def browseFiles_reference():
    # Change the default to "/" and allow user input. 
    filename = filedialog.askopenfilename(initialdir = "/Users/jenniliu/Desktop/BIOINF576/TCGA_Matchmaker/TCGA_code/example_test")
      
    # Display the file path
    ent_ref_profile.insert(0, filename)


def browseFiles_TCGA():
    # Change the default to "/" and allow user input. 
    filename = filedialog.askopenfilename(initialdir = "/Users/jenniliu/Desktop/BIOINF576/TCGA_Matchmaker/TCGA_code/example_test")
      
    # Display the file path
    ent_tcga_profile.insert(0, filename)


def analysis_test():
    # Pre-process data:
    #input_profile = m.read_expr_profile(input_profile_path)
    path = "/Users/jenniliu/Desktop/BIOINF576/TCGA_Matchmaker/TCGA_code/example_test/example_input_sample.csv"
    input_profile = pd.read_csv(path, sep = ";", encoding="UTF-8")
    # Average the duplicate values
    #input_profile = input_profile.groupby('symbol', as_index=False).mean()


    lbl_output_window["text"] = str(input_ref_profile_var.get())
 #'HI' +str(missing_values_var.get())



def analysis(add_missing = False, output = True):#(input_profile_path, input_sample_data_path, add_missing = False, output = True):
    #print('hello')
    
    # IF INPUT PATHS ARE EMPTY RETURN MESSAGE TO THE OUTPUT WINDOW!


    # Pre-process data:
    #input_profile = m.read_expr_profile(input_profile_path)
    path = input_ref_profile_var.get()
    input_profile = pd.read_csv(path, sep = ";", encoding="UTF-8")
    # Average the duplicate values
    input_profile = input_profile.groupby('symbol', as_index=False).mean()

    #input_sample_data = m.read_TCGA_sample(input_sample_data_path)
    path = input_tcga_profile_var.get()
    input_sample_data = pd.read_csv(path, sep = ";", encoding="UTF-8")
    # Average the duplicate values
    input_sample_data = input_sample_data.groupby('symbol', as_index=False).mean()

    profile_level = input_profile.iloc[:,1]
    is_nan = np.all(np.isnan(profile_level))
    print("IS THERE is_nan?", is_nan)

    # Clean data:
    profile, sample, missing_TCGA = m.check_TCGA(input_profile, input_sample_data, add_missing, output)
    input_ref_profile_df, input_tcga_profile_df, missing_reference = m.check_profile(profile, sample, add_missing, output)

    # "profile" and "sample" are now pre-processed and clean.


    # MAYBE: YOU SHOULD CHECK FOR ISNAN OR ISINF HERE FOR input_ref_profile_df, input_tcga_profile_df

    # Compute distance:
    distance = m.compute_distance(input_ref_profile_df, input_tcga_profile_df)

    #print("The input expression levels show a correlation value of:", distance, "when zero expression levels are added.")
    print("Genes that are missing in the reference profile are:\n", missing_TCGA)
    print("Genes that are missing in the TCGA profile are:\n", missing_reference)

    #lbl_output_window["text"] = "HELLO"
    lbl_output_window["text"] = f"The input expression levels show a correlation value of: {distance} when zero expression levels are added."




label = tk.Label(text="TCGA Matchmaker")
label.pack()

############################

# First Frame: Load input files
# command=partial(read_expr_profile, filename)
frm_load_input = tk.Frame()
frm_load_input.pack(fill=tk.X, ipadx=5, ipady=5)

btn_input_ref = tk.Button(master=frm_load_input, text="Input Reference Profile", command=browseFiles_reference)
btn_input_ref.grid(row=1, column=0, sticky='e')
btn_input_TCGA = tk.Button(master=frm_load_input, text="Input TCGA Profile", command=browseFiles_TCGA)
btn_input_TCGA.grid(row=1, column=5, sticky='e')

# Display and save the file paths.
frm_input_files = tk.Frame(relief=tk.SUNKEN, borderwidth=3)
frm_input_files.pack()

lbl_ref_profile = tk.Label(master=frm_input_files, text="Input Reference Profile:")
ent_ref_profile = tk.Entry(master=frm_input_files, width=50, textvariable = input_ref_profile_var)
lbl_ref_profile.grid(row=0, column=0, sticky="e")
ent_ref_profile.grid(row=0, column=1)

lbl_tcga_profile = tk.Label(master=frm_input_files, text="Input TCGA Profile:")
ent_tcga_profile = tk.Entry(master=frm_input_files, width=50, textvariable = input_tcga_profile_var)
lbl_tcga_profile.grid(row=1, column=0, sticky="e")
ent_tcga_profile.grid(row=1, column=1)

btn_enter_files = tk.Button(master=frm_input_files, text = "Enter", command=enter_profiles)
btn_enter_files.grid(row=2, column=2)


label.pack()
# Set Parameters Frame
label.pack()
label = tk.Label(text="Enter Parameters")
label.pack()

frm_param = tk.Frame(relief=tk.SUNKEN, borderwidth=3)
frm_param.pack()

lbl_missing_values = tk.Label(master=frm_param, text="Missing Values:")
ent_missing_values = tk.Entry(master=frm_param, width=50, textvariable = missing_values_var)
lbl_missing_values.grid(row=0, column=0, sticky="e")
ent_missing_values.grid(row=0, column=1)


lbl_distance = tk.Label(master=frm_param, text="Distance Metrics:")
ent_distance = tk.Entry(master=frm_param, width=50, textvariable = distance_metrics_var)
lbl_distance.grid(row=1, column=0, sticky="e")
ent_distance.grid(row=1, column=1)

lbl_output = tk.Label(master=frm_param, text="Show Output:")
ent_output = tk.Entry(master=frm_param, width=50, textvariable = show_output_var)
lbl_output.grid(row=2, column=0, sticky="e")
ent_output.grid(row=2, column=1)

# Enter all parameters at once
btn_enter_parameters = tk.Button(master=frm_param, text="Enter", command=enter_entry)
btn_enter_parameters.grid(row=3, column=2)



# Analysis Button Frame
frm_Analysis = tk.Frame()
frm_Analysis.pack(fill=tk.X, ipadx=5, ipady=5)

# btn_Analysis = tk.Button(master=frm_Analysis, text="Start Analysis", command=partial(analysis, 
#                                                                                      input_ref_profile_var.get(), 
#                                                                                      input_tcga_profile_var.get(), 
#                                                                                      missing_values_var.get(),
#                                                                                      show_output_var.get()))

btn_Analysis = tk.Button(master=frm_Analysis, text="Start Analysis", command=partial(analysis, missing_values_var.get(),
                                                                                      show_output_var.get()))

#btn_Analysis = tk.Button(master=frm_Analysis, text="Start Analysis", command=analysis_test)

btn_Analysis.pack(side=tk.LEFT, ipadx=10)


# Output Window Frame
label.pack()
label = tk.Label(text="Output Window")
label.pack()

lbl_output_window = tk.Label(
    text="",
    fg="black",
    bg="white",
    relief=tk.SUNKEN, 
    borderwidth=3
)
lbl_output_window.pack(fill=tk.X, ipadx=50, ipady=50)















window.mainloop()
