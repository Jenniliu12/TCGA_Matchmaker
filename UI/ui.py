import tkinter as tk
from tkinter import ttk
from tkinter import *
from tkinter import filedialog
import pandas as pd
import numpy as np
from functools import partial
import csv
import match_computation as m


# START WINDOW EVENT LOOP:
window = tk.Tk()
#window.geometry("750x750")



input_ref_profile_var = tk.StringVar()
input_tcga_profile_var = tk.StringVar()
output_path_var = tk.StringVar()
missing_values_var = tk.StringVar()
method_var = tk.StringVar()
show_output_var = tk.StringVar()
threshold_var = tk.StringVar()
output_text1 = tk.StringVar()
output_text2 = tk.StringVar()


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

def browseFiles_output():
    # Change the default to "/" and allow user input. 
    filename = filedialog.askdirectory(initialdir = "/Users/jenniliu/Desktop/BIOINF576/TCGA_Matchmaker/TCGA_code/example_test")
      
    # Display the file path
    ent_output_path.insert(0, filename)


def analysis():
    
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

    # Clean data:
    profile, sample, missing_TCGA = m.check_TCGA(input_profile, input_sample_data, missing_values_var.get() == "True", show_output_var.get() == "True")
    profile, sample, missing_reference = m.check_profile(profile, sample, missing_values_var.get() == "True", show_output_var.get() == "True")

    # Normalize data sets
    profile = m.normalize_profile(profile, method_var.get())
    sample = m.normalize_profile(sample, method_var.get())

    # "profile" and "sample" are now pre-processed and clean.

    # MAYBE: YOU SHOULD CHECK FOR ISNAN OR ISINF HERE FOR profile, input_tcga_profile_df

    # ACTUAL ANALYSIS STEP:
    # (1) Compute distance:
    distance = m.compute_distance(profile, sample)
    file_path = output_path_var.get() + "/OUTPUT_TEST.csv"
    with open(file_path, 'w+', newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ' ')
        text = ""
        text += "Missing Genes reference profile are:\n" #+ missing_TCGA
        for gene in missing_reference:
            text += f"\n{gene}"

        text += "\nMissing Genes TCGA profile are:" #+ missing_reference
        for gene in missing_TCGA:
            text += f"\n{gene}"

        my_writer.writerow("Hellloo")
        my_writer.writerow(text)

        # TEST DOES IT WORK AT ALL?
        print("HELLLLLLLOOOOO")

    if show_output_var.get() == "True":
        print(text)

    # Parameter and setting output:
    if show_output_var.get() == "True":
        text = ""
        text += f"Function output is shown in the Terminal."
    else:
        text = ""
        text += f"Function output is not shown in the Terminal."

    if missing_values_var.get() == "True":
        text += f"\nMissing genes were added with zero expression values. No genes were dropped from the dataset."
    else:
        text += f"\nMissing genes were dropped from the dataset."

    text += f"\nThe threshold sensitivity that determines similiarity between gene expression levels is {threshold_var.get()}."
    text += f"\nThe chosen normalization method is {method_var.get()}."
    text += f"\nIf any genes are missing, they are collected in an external file."
    output_text1.set(text)

    # Analysis output:
    text = ""
    text += f"The input expression levels have a pearson correlation coefficient of: {distance}."
    text += f"\nGenes with similiar expression levels threshold are collected in an external file."
    output_text2.set(text)

    # (2) Bar Chart:
    gene_ratio, sim_genes = m.expression_analysis(profile, sample)
    m.gene_bar_chart(gene_ratio)



title = tk.Label(text="TCGA Matchmaker", font=('Arial', 25))
title.pack()
descript = tk.Label(text="This Tool compares Gene Expression Profiles to find cancerous Biomarkers.")
descript.pack()
descript2 = tk.Label(text="Follow the following instructions to run the tool:")
descript2.pack()
empty = tk.Label()
empty.pack()
instructions = tk.Label(text="1: Enter input files (.csv) that are formatted: 2 Columns ('symbol', 'value'). Enter output file path.")
instructions.pack(pady=0, side= TOP, anchor="w")
instructions2 = tk.Label(text="2: Enter parameter values and settings.")
instructions2.pack(pady=0, side= TOP, anchor="w")
instructions3 = tk.Label(text="3: Click 'Start Analysis' to run the tool.")
instructions3.pack(pady=0, side= TOP, anchor="w")
empty = tk.Label()
empty.pack()
############################

# First Frame: Load input files
frm_load_input = tk.Frame()
frm_load_input.pack(fill=tk.X, ipadx=5, ipady=5)

btn_input_ref = tk.Button(master=frm_load_input, text="Input Reference Profile", command=browseFiles_reference)
btn_input_ref.grid(row=1, column=0, sticky='e')
btn_input_TCGA = tk.Button(master=frm_load_input, text="Input TCGA Profile", command=browseFiles_TCGA)
btn_input_TCGA.grid(row=1, column=5, sticky='e')

# Enter output path
btn_input_ref = tk.Button(master=frm_load_input, text="Choose Output Path", command=browseFiles_output)
btn_input_ref.grid(row=1, column=10, sticky='e')

# Display and save the file paths.
frm_input_files = tk.Frame(relief=tk.SUNKEN, borderwidth=3)
frm_input_files.pack()

lbl_ref_profile = tk.Label(master=frm_input_files, text="Input Reference Profile:")
ent_ref_profile = tk.Entry(master=frm_input_files, width=70, textvariable = input_ref_profile_var)
lbl_ref_profile.grid(row=0, column=0, sticky="e")
ent_ref_profile.grid(row=0, column=1)

lbl_tcga_profile = tk.Label(master=frm_input_files, text="Input TCGA Profile:")
ent_tcga_profile = tk.Entry(master=frm_input_files, width=70, textvariable = input_tcga_profile_var)
lbl_tcga_profile.grid(row=1, column=0, sticky="e")
ent_tcga_profile.grid(row=1, column=1)

# Display and save the chosen output file path.
lbl_ref_profile = tk.Label(master=frm_input_files, text="Output Path:")
ent_output_path = tk.Entry(master=frm_input_files, width=70, textvariable = output_path_var)
lbl_ref_profile.grid(row=2, column=0, sticky="e")
ent_output_path.grid(row=2, column=1)


title.pack()
# Set Parameters Frame
title.pack()
label = tk.Label(text="Enter Parameters")
label.pack()

frm_param = tk.Frame(relief=tk.SUNKEN, borderwidth=3)
frm_param.pack()

lbl_missing_values = tk.Label(master=frm_param, text="Missing Values:")
ent_missing_values = tk.Entry(master=frm_param, width=30, textvariable = missing_values_var)
lbl_missing_values.grid(row=0, column=0, sticky="e")
ent_missing_values.grid(row=0, column=1)

lbl_output = tk.Label(master=frm_param, text="Show Output:")
ent_output = tk.Entry(master=frm_param, width=30, textvariable = show_output_var)
lbl_output.grid(row=1, column=0, sticky="e")
ent_output.grid(row=1, column=1)

lbl_distance = tk.Label(master=frm_param, text="Normalization Method:")
ent_distance = tk.Entry(master=frm_param, width=30, textvariable = method_var)
lbl_distance.grid(row=2, column=0, sticky="e")
ent_distance.grid(row=2, column=1)

lbl_distance = tk.Label(master=frm_param, text="Threshold value (0-1):")
ent_distance = tk.Entry(master=frm_param, width=30, textvariable = threshold_var)
lbl_distance.grid(row=3, column=0, sticky="e")
ent_distance.grid(row=3, column=1)


# Analysis Button Frame
frm_Analysis = tk.Frame()
frm_Analysis.pack(fill=tk.X, ipadx=5, ipady=5)
btn_Analysis = tk.Button(master=frm_Analysis, text="Start Analysis", command=analysis)
btn_Analysis.pack(side=tk.LEFT, ipadx=10)



# Output Window Frame
label = tk.Label(text="Summary: Configurations and Parameters")
label.pack()

lbl_output_window1 = tk.Label(
    textvariable=output_text1,
    fg="black",
    bg="white",
    relief=tk.SUNKEN, 
    borderwidth=1,
    #anchor='w', 
    #justify='left'
)
lbl_output_window1.pack(fill=tk.X,ipadx=30, ipady=30, side=TOP, anchor="w")

label = tk.Label(text="Analysis Output")
label.pack()

lbl_output_window2 = tk.Label(
    textvariable=output_text2,
    fg="black",
    bg="white",
    relief=tk.SUNKEN, 
    borderwidth=3
)
lbl_output_window2.pack(fill=tk.X, ipadx=30, ipady=30)











window.mainloop()
