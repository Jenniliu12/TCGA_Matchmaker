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
window.geometry("750x750")
window.config(background="#C7DEFF")


input_ref_profile_var = tk.StringVar()
input_tcga_profile_var = tk.StringVar()
missing_values_var = tk.StringVar()
method_var = tk.StringVar()
show_output_var = tk.StringVar()
info_text = tk.StringVar()
output_text = tk.StringVar()
window_text = ""

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


    #lbl_output_window["text"] = input_ref_profile_var.get()
    
    #'HI' +str(missing_values_var.get())



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

    
    gene_ratio, sim_genes = m.expression_analysis(profile, sample)
    m.gene_bar_chart(gene_ratio)

    # Compute distance:
    distance = m.compute_distance(profile, sample)

    if show_output_var.get() == "True":
        text = ""
        text += "Genes that are missing in the reference profile are:\n" #+ missing_TCGA
        for gene in missing_reference:
            text += f"\n{gene}"

        text += "Genes that are missing in the TCGA profile are:\n" #+ missing_reference
        for gene in missing_TCGA:
            text += f"\n{gene}"

        info_text.set(text)

    text = ""
    text += f"\nThe input expression levels show a correlation value of: {distance}."
    text += f"\nMissing expression values were added: {missing_values_var.get()}."
    text += f"\nShow output?: {show_output_var.get()}."
    output_text.set(text)




label = tk.Label(text="TCGA Matchmaker", background="#C7DEFF", font=('Arial', 25))
label.pack()
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
ent_ref_profile = tk.Entry(master=frm_input_files, width=70, textvariable = input_ref_profile_var)
lbl_ref_profile.grid(row=0, column=0, sticky="e")
ent_ref_profile.grid(row=0, column=1)

lbl_tcga_profile = tk.Label(master=frm_input_files, text="Input TCGA Profile:")
ent_tcga_profile = tk.Entry(master=frm_input_files, width=70, textvariable = input_tcga_profile_var)
lbl_tcga_profile.grid(row=1, column=0, sticky="e")
ent_tcga_profile.grid(row=1, column=1)


label.pack()
# Set Parameters Frame
label.pack()
label = tk.Label(text="Enter Parameters", background="#C7DEFF")
label.pack()

frm_param = tk.Frame(relief=tk.SUNKEN, borderwidth=3)
frm_param.pack()

lbl_missing_values = tk.Label(master=frm_param, text="Add Missing Values:")
ent_missing_values = tk.Entry(master=frm_param, width=70, textvariable = missing_values_var)
lbl_missing_values.grid(row=0, column=0, sticky="e")
ent_missing_values.grid(row=0, column=1)

lbl_output = tk.Label(master=frm_param, text="Show Output:")
ent_output = tk.Entry(master=frm_param, width=70, textvariable = show_output_var)
lbl_output.grid(row=1, column=0, sticky="e")
ent_output.grid(row=1, column=1)

lbl_distance = tk.Label(master=frm_param, text="Normalization Method:")
ent_distance = tk.Entry(master=frm_param, width=70, textvariable = method_var)
lbl_distance.grid(row=2, column=0, sticky="e")
ent_distance.grid(row=2, column=1)



# Analysis Button Frame
frm_Analysis = tk.Frame()
frm_Analysis.pack(fill=tk.X, ipadx=5, ipady=5)


btn_Analysis = tk.Button(master=frm_Analysis, text="Start Analysis", command=partial(analysis))

btn_Analysis.pack(side=tk.LEFT, ipadx=10)





# Output Window Frame
label.pack()
label = tk.Label(text="Output Window", background="#C7DEFF")
label.pack()


# Scrollbar Widget


# entr_output_window = tk.Entry(
#     textvariable=info_text,
#     fg="black",
#     bg="white",
#     relief=tk.SUNKEN, 
#     borderwidth=3,
# )
# v = Scrollbar(entr_output_window, orient='vertical')
# v2 = Scrollbar(entr_output_window, orient='horizontal')
# v.pack(side=RIGHT, fill='y')
# v2.pack(side=BOTTOM, fill='y')

#entr_output_window.pack(fill=tk.X, ipadx=100, ipady=100)


lbl_output_window2 = tk.Label(
    textvariable=output_text,
    fg="black",
    bg="white",
    relief=tk.SUNKEN, 
    borderwidth=3
)
lbl_output_window2.pack(fill=tk.X, ipadx=50, ipady=50)














window.mainloop()
