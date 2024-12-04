import tkinter as tk
from tkinter import *
# AA weights https://web.expasy.org/findmod/findmod_masses.html#AA
complement_dictionary = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N", "X":"X"}
dna_to_protein = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

amino_acid_weights = {
    'A': 71.0788,
    'R': 156.1875,
    'N': 114.1038,
    'D': 115.0886,
    'C': 103.1388,
    'E': 129.1155,
    'Q': 128.1307,
    'G': 57.0519,
    'H': 137.1411,
    'I': 113.1594,
    'L': 113.1594,
    'K': 128.1741,
    'M': 131.1926,
    'F': 147.1766,
    'P': 97.1167,
    'S': 87.0782,
    'T': 101.1051,
    'W': 186.2132,
    'Y': 163.176,
    'V': 99.1326
}




class MultiPageApp(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("DNA and Protein property calculator")
        self.geometry("500x600")

        # Container to hold all frames
        container = tk.Frame(self)
        container.pack(fill="both", expand=True)

        self.frames = {}

        # Initialize all frames and add them to the container
        for F in (StartPage, PageOne, Pagetwo, Pagethree, Pagefour):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    def show_frame(self, page):
        frame = self.frames[page]
        frame.tkraise()


class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)

        label = tk.Label(self, text="Choose what do you want to do?", font=("Helvetica", 14))
        label.grid(row = 0, column = 0, sticky = "W", pady = 2, columnspan = 3)

        next1_button = tk.Button(self, text="DNA sequence operations", bg = "light goldenrod", width = 20, command=lambda: controller.show_frame(PageOne))
        next1_button.grid(row = 1, column = 0, sticky = "W", pady = 5, padx = 5)
        next2_button = tk.Button(self, text="protein A280 to conc.", bg="yellow", width = 20, command=lambda: controller.show_frame(Pagetwo))
        next2_button.grid(row = 1, column = 1, sticky = "W", pady = 5, padx = 5)

        next3_button = tk.Button(self, text="conc mg/ml to uM", bg="NavajoWhite2", width = 20, command=lambda: controller.show_frame(Pagethree))
        next3_button.grid(row = 2, column = 0, sticky = "W", pady = 5, padx = 5)

        next4_button = tk.Button(self, text="conc uM to mg/ml", bg="gold", width = 20, command=lambda: controller.show_frame(Pagefour))
        next4_button.grid(row = 2, column = 1, sticky = "W", pady = 5, padx = 5)


class PageOne(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)

        label1 = tk.Label(self, text="Find complement, reverse complement", font=("Helvetica", 12))
        label1.grid(row = 0, column = 0, sticky = "W", pady = 2, columnspan = 3)
        label2 = tk.Label(self, text="or translate your DNA", font=("Helvetica", 12))
        label2.grid(row=1, column=0, sticky="W", pady=2, padx=2, columnspan=3)
        label3 = tk.Label(self, text="Paste your sequence in the window below", font=("Helvetica", 12))
        label3.grid(row=2, column=0, sticky="W", pady=2, padx=2, columnspan=3)

        seq_entry = tk.Text(self, font=('calibre', 10, 'normal'), width = 50, height = 8)
        seq_entry.grid(row = 3, column = 0, sticky = "W", pady = 2, padx = 10, columnspan = 5)
        label4 = tk.Label(self, text="Choose what do you want to do with it", font=("Helvetica", 12))
        label4.grid(row = 4, column = 0, sticky = "W", pady = 2, columnspan = 4)
        complement_button = tk.Button(self, text="Complement", bg = "NavajoWhite2", command= lambda: getting_value_complement())
        complement_button.grid(row = 5, column = 0, sticky = "W", pady = 2, padx = 5)
        reverse_complement_button = tk.Button(self, text="Reverse complement", bg="NavajoWhite1", command= lambda: getting_value_reverse_complement())
        reverse_complement_button.grid(row=5, column=1, sticky="W", pady=2, padx=2)
        translate_button = tk.Button(self, text="Translate", bg="LightGoldenrod1", command= lambda: translate_seq())
        translate_button.grid(row=5, column=2, sticky="W", pady=2, padx=2)
        clear_button = tk.Button(self, text="Clear", bg="wheat1", command = lambda:clearall())
        clear_button.grid(row=5, column=3, sticky="W", pady=2, padx=2)
        label5 = tk.Label(self, text="Output:", font=("Helvetica", 12))
        label5.grid(row = 6, column = 0, sticky = "W", pady = 2, padx = 5, columnspan = 2)

        seq_out = tk.Text(self, font=('calibre', 10, 'normal'), width=50, height=8)
        seq_out.grid(row=7, column=0, sticky="W", pady=2, padx=10, columnspan=5)
        back_button = tk.Button(self, text="Back", width = 8, bg = "wheat1", command=lambda: controller.show_frame(StartPage))
        back_button.grid(row = 8, column = 2, sticky = "E", pady = 5, padx = 8, columnspan = 3)

        def clearall():
            seq_out.delete('1.0', END)
            seq_entry.delete('1.0', END)

        def getting_value_complement():
            sequence = str(seq_entry.get(1.0, "end-1c"))
            sequence1 = wrong_element_remover(sequence.upper())
            new_seq = ""
            for a in sequence1:
                new_seq += complement_dictionary[a]
            seq_out.delete('1.0', END)
            seq_out.insert(END, new_seq)

        def getting_value_reverse_complement():
            sequence = str(seq_entry.get(1.0, "end-1c"))
            sequence2 = wrong_element_remover(sequence.upper())
            sequence1 = sequence2[::-1]
            new_seq = ""
            for a in sequence1:
                new_seq += complement_dictionary[a]
            seq_out.delete('1.0', END)
            seq_out.insert(END, new_seq)

        def translate_seq():
            text_output = """"""
            seq_list = []
            label_list = ["frame 1", "frame 2", "frame 3", "frame 1 reverse", "frame 2 reverse", "frame 3 reverse"]
            sequence = str(seq_entry.get(1.0, "end-1c"))
            sequence_upcase = sequence.upper()
            sequence_frame1 = wrong_element_remover(sequence_upcase)
            sequence_frame2 = sequence_frame1[1:]
            sequence_frame3 = sequence_frame1[2:]
            sequence_frame1_inv = sequence_frame1[::-1]
            reverse_comp_frw1 = "".join([complement_dictionary[a] for a in list(sequence_frame1_inv)])
            reverse_comp_frw2 = reverse_comp_frw1[1:]
            reverse_comp_frw3 = reverse_comp_frw1[2:]
            seq_list.extend([sequence_frame1, sequence_frame2, sequence_frame3, reverse_comp_frw1, reverse_comp_frw2, reverse_comp_frw3])
            for elem, label in zip(seq_list, label_list):
                translated_elem = translator(elem)
                text_output += f"{label}: {translated_elem}\n"
            seq_out.delete('1.0', END)
            seq_out.insert(END, text_output)

        def translator(sequence):
            triplets_preliminary = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
            triplets = [a for a in triplets_preliminary if len(a) == 3]
            prot_seq = "".join([dna_to_protein[triplet] if triplet in dna_to_protein.keys() else "X" for triplet in triplets])
            return prot_seq

        def wrong_element_remover(sequence):
            sequence1 = sequence.replace("\n", "")
            sequence2 = "".join([a for a in list(sequence1) if a in complement_dictionary.keys()])
            return sequence2





class Pagetwo(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)

        label1 = tk.Label(self, text="Calculate concentration from A280", font=("Helvetica", 14))
        label1.grid(row = 0, column = 0, sticky = "W", pady = 2, padx = 5, columnspan = 6)
        label2 = tk.Label(self, text="Input below either protein sequence or extinction coefficient.", font=("Helvetica", 11))
        label2.grid(row=1, column=0, sticky="W", pady=2, padx = 5, columnspan=6)
        label3 = tk.Label(self, text="You can leave molecular weight box blank", font=("Helvetica", 11))
        label3.grid(row=2, column=0, sticky="W", pady=2, padx = 5, columnspan=6)
        label4 = tk.Label(self, text="If you want concentration in mg/ml as well, you need to", font=("Helvetica", 11))
        label4.grid(row=3, column=0, sticky="W", pady=2, padx=5, columnspan=6)
        label5 = tk.Label(self, text="provide molecular weight of the protein in Da (g/mol),", font=("Helvetica", 11))
        label5.grid(row=4, column=0, sticky="W", pady=2, padx=5, columnspan=6)
        label6 = tk.Label(self, text="or just use the protein sequence option", font=("Helvetica", 11))
        label6.grid(row=5, column=0, sticky="W", pady=2, padx=5, columnspan=6)
        label7 = tk.Label(self, text="", font=("Helvetica", 5))
        label7.grid(row=6, column=0, sticky="W", pady=2, padx = 5, columnspan=6)
        label8 = tk.Label(self, text = "Type in your A280 value in the box below", font=("Helvetica", 12))
        label8.grid(row=7, column=0, sticky="W", pady=2, padx = 5, columnspan=6)
        a280entry = tk.Entry(self, font = ('calibre',10,'normal'), width = 10, bd = 5)
        a280entry.grid(row=8, column=0, sticky="W", pady=2, padx = 5, columnspan=6)
        label9 = tk.Label(self, text="Input the extinction coefficient (e) and (optional) molecular weight in Da", font=("Helvetica", 11))
        label9.grid(row=9, column=0, sticky="W", pady=2, padx=5, columnspan=6)
        label_ext = tk.Label(self, text="e:", font=("Helvetica", 11))
        label_ext.grid(row=10, column=0, sticky="W", pady=2, padx=5, columnspan=7)
        ext_coef = tk.Entry(self, font = ('calibre',10,'normal'), width = 10)
        ext_coef.grid(row=10, column=1, sticky="W", pady=2, padx= 2, columnspan=7)
        label_mwt = tk.Label(self, text="Mwt (Da):", font=("Helvetica", 11))
        label_mwt.grid(row=10, column=2, sticky="W", pady=2, padx=10, columnspan=7)
        mwt_entry = tk.Entry(self, font=('calibre', 10, 'normal'), width=10)
        mwt_entry.grid(row=10, column=3, sticky="W", pady=2, padx=2, columnspan=7)
        label9 = tk.Label(self, text="input the protein sequence below", font=("Helvetica", 11))
        label9.grid(row=11, column=0, sticky="W", pady=2, padx=5, columnspan=7)
        seq_entry = tk.Text(self, font=('calibre', 10, 'normal'), width=50, height=5)
        seq_entry.grid(row=12, column=0, sticky="W", pady=2, padx=10, columnspan=7)
        label_output = tk.Label(self, text="Output:", font=("Helvetica", 11))
        label_output.grid(row=13, column=0, sticky="W", pady=2, padx=10, columnspan=7)
        out_box = tk.Text(self, font=('calibre', 10, 'normal'), width=50, height=5)
        out_box.grid(row=14, column=0, sticky="W", pady=2, padx=10, columnspan=7)
        calculate_button = tk.Button(self, text="Calculate", bg="bisque2", width=10, command=lambda : get_values_a280())
        calculate_button.grid(row=15, column=0, sticky="W", pady=2, padx = 10,  columnspan=6)
        clear_button = tk.Button(self, text="Clear", bg="light cyan", width=8, command= lambda : clearall_1())
        clear_button.grid(row=15, column=2, sticky="W", pady=2, padx = 1, columnspan=7)
        back_button = tk.Button(self, text="Back", bg="NavajoWhite2", width=8, command=lambda: controller.show_frame(StartPage))
        back_button.grid(row=15, column=3, sticky="W", pady=2, columnspan=6, padx=2)

        def get_values_a280():
            A280 = float(a280entry.get())
            sequence3 = str(seq_entry.get(1.0, "end-1c"))
            sequence = sequence3.replace("\n", "")

            extinct_coef_value = ext_coef.get()
            mwt_value = mwt_entry.get()
            output_text = ""
            if extinct_coef_value != "":
                conc_uM = A280/(float(extinct_coef_value))*1000000
                output_text += f"molar concentration is: \n {round(conc_uM, 2)} uM \n"
                if mwt_value != "":
                    conc_mg_ml = conc_uM*float(mwt_value)/1000000
                    output_text += f"Concentration mg/ml: {round(conc_mg_ml, 2)} \n"
            elif extinct_coef_value == "":
                mwt_calculated = 18.02
                sequence1 = "".join([a for a in list(sequence) if a in amino_acid_weights.keys()])
                print(sequence1)
                extinction_calculated = sequence1.count("W") * 5500 + sequence1.count("Y") * 1490 + (sequence1.count("C") // 2) * 125
                extinction_calculated1 = sequence1.count("W") * 5500 + sequence1.count("Y") * 1490

                for aa in sequence1:
                    mwt_calculated += amino_acid_weights[aa]
                conc_uM = A280/extinction_calculated*1000000
                conc_mg_ml = conc_uM*mwt_calculated/1000000
                output_text += f" concentration uM {round(conc_uM, 2)}, \n concentration mg/ml {round(conc_mg_ml, 2)}, \n Molecular weight {round(mwt_calculated, 2)} Da, \n Extinction coefficient if cysteines form cystines {round(extinction_calculated, 2)} \n Extinction coefficient if cysteines arr reduced {round(extinction_calculated1, 2)}"
            out_box.delete('1.0', END)
            out_box.insert(END, output_text)

        def clearall_1():
            out_box.delete('1.0', END)
            a280entry.delete('0', END)
            ext_coef.delete('0', END)
            mwt_entry.delete('0', END)
            seq_entry.delete('1.0', END)


class Pagethree(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)

        label1 = tk.Label(self, text="mg/ml to uM converter", font=("Helvetica", 16))
        label1.grid(row = 0, column = 0, sticky = "W", pady = 2, padx = 5, columnspan = 6)
        label2 = tk.Label(self, text="Input concentration in mg/ml in the box below", font=("Helvetica", 12))
        label2.grid(row = 1, column = 0, sticky = "W", pady = 2, padx = 5, columnspan = 6)
        mg_ml_conc_entr = tk.Entry(self, font = ('calibre',10,'normal'), width = 10, bd = 5)
        mg_ml_conc_entr.grid(row = 2, column = 0, sticky = "W", pady = 2, padx = 10, columnspan = 6)
        label3 = tk.Label(self, text="Input either molecular weight or protein sequence in the \n corresponding box below, choose the right box", font=("Helvetica", 12))
        label3.grid(row = 3, column = 0, sticky = "W", pady = 2, padx = 5, columnspan = 6)
        label4 = tk.Label(self, text="Input molecular weight in Da below:", font=("Helvetica", 12))
        label4.grid(row = 4, column = 0, sticky = "W", pady = 2, padx = 5, columnspan = 6)
        mwt_box_entr = tk.Entry(self, font=('calibre', 10, 'normal'), width=15)
        mwt_box_entr.grid(row=5, column=0, sticky="W", pady=2, padx=10, columnspan=6)

        label5 = tk.Label(self, text="Input protein sequence below:", font=("Helvetica", 12))
        label5.grid(row=6, column=0, sticky="W", pady=2, padx=5, columnspan=6)
        seq_protentry = tk.Text(self, font=('calibre', 10, 'normal'), width=50, height=8)
        seq_protentry.grid(row=7, column=0, sticky="W", pady=2, padx=10, columnspan=5)
        label6 = tk.Label(self, text="Output:", font=("Helvetica", 12))
        label6.grid(row=8, column=0, sticky="W", pady=2, padx=5, columnspan=6)
        seq_protout = tk.Text(self, font=('calibre', 10, 'normal'), width=50, height=5)
        seq_protout.grid(row=9, column=0, sticky="W", pady=2, padx=10, columnspan=5)
        back_button = tk.Button(self, text="Back", width = 10, bg = "bisque2", command=lambda: controller.show_frame(StartPage))
        back_button.grid(row = 10, column = 3, sticky = "W", pady = 2, ipadx = 10, columnspan = 6)
        calculate_button1 = tk.Button(self, text="Calculate", bg="LemonChiffon3", width=10, command=lambda:get_out_entr())
        calculate_button1.grid(row=10, column=0, sticky="W", pady=2, padx=10, columnspan=6)
        clear_button = tk.Button(self, text="Clear", bg="light cyan", width=8, command= lambda : clearall_2())
        clear_button.grid(row=10, column=1, sticky="NE", pady=2, padx=5, columnspan=2)

        def get_out_entr():
            out_text_mg = ""
            sequence4 = str(seq_protentry.get(1.0, "end-1c"))
            sequence5 = sequence4.replace("\n", "")
            conc_mg_ml_val = float(mg_ml_conc_entr.get())
            if sequence5 != "":
                sequence6 = "".join([a for a in list(sequence5) if a in amino_acid_weights.keys()])
                mwt_calculated = 18.02 + sum([amino_acid_weights[a] for a in sequence6])
                res = conc_mg_ml_val/float(mwt_calculated)*1000000
                out_text_mg += f"Concentration is {round(res, 2)} uM, \nMolecular weight is {round(mwt_calculated, 2)} Da"
            elif sequence5 == "":
                mwt_provided = float(mwt_box_entr.get())
                res1 = conc_mg_ml_val / mwt_provided * 1000000
                out_text_mg += f"Concentration is {round(res1, 2)} uM"
            seq_protout.delete('1.0', END)
            seq_protout.insert(END, out_text_mg)

        def clearall_2():
            seq_protout.delete('1.0', END)
            seq_protentry.delete('1.0', END)
            mwt_box_entr.delete('0', END)
            mg_ml_conc_entr.delete('0', END)






class Pagefour(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)

        label1 = tk.Label(self, text="uM to mg/ml converter", font=("Helvetica", 16))
        label1.grid(row = 0, column = 0, sticky = "W", pady = 2, padx = 5, columnspan = 6)
        label2 = tk.Label(self, text="Input concentration in uM in the box below", font=("Helvetica", 12))
        label2.grid(row = 1, column = 0, sticky = "W", pady = 2, padx = 5, columnspan = 6)
        uM_conc_entr = tk.Entry(self, font = ('calibre',10,'normal'), width = 10, bd = 5)
        uM_conc_entr.grid(row = 2, column = 0, sticky = "W", pady = 2, padx = 10, columnspan = 6)
        label3 = tk.Label(self, text="Input either molecular weight or protein sequence in the \n corresponding box below, choose the right box", font=("Helvetica", 12))
        label3.grid(row = 3, column = 0, sticky = "W", pady = 2, padx = 5, columnspan = 6)
        label4 = tk.Label(self, text="Input molecular weight in Da below:", font=("Helvetica", 12))
        label4.grid(row = 4, column = 0, sticky = "W", pady = 2, padx = 5, columnspan = 6)
        mwt_box_entr = tk.Entry(self, font=('calibre', 10, 'normal'), width=15)
        mwt_box_entr.grid(row=5, column=0, sticky="W", pady=2, padx=10, columnspan=6)

        label5 = tk.Label(self, text="Input protein sequence below:", font=("Helvetica", 12))
        label5.grid(row=6, column=0, sticky="W", pady=2, padx=5, columnspan=6)
        seq_protentry = tk.Text(self, font=('calibre', 10, 'normal'), width=50, height=8)
        seq_protentry.grid(row=7, column=0, sticky="W", pady=2, padx=10, columnspan=5)
        label6 = tk.Label(self, text="Output:", font=("Helvetica", 12))
        label6.grid(row=8, column=0, sticky="W", pady=2, padx=5, columnspan=6)
        seq_protout = tk.Text(self, font=('calibre', 10, 'normal'), width=50, height=5)
        seq_protout.grid(row=9, column=0, sticky="W", pady=2, padx=10, columnspan=5)
        back_button = tk.Button(self, text="Back", width = 10, bg = "bisque2", command=lambda: controller.show_frame(StartPage))
        back_button.grid(row = 10, column = 3, sticky = "W", pady = 2, ipadx = 10, columnspan = 6)
        calculate_button1 = tk.Button(self, text="Calculate", bg="thistle1", width=10, command=lambda:get_out_entr())
        calculate_button1.grid(row=10, column=0, sticky="W", pady=2, padx=10, columnspan=6)
        clear_button = tk.Button(self, text="Clear", bg="light cyan", width=8, command= lambda : clearall_3())
        clear_button.grid(row=10, column=1, sticky="NE", pady=2, padx=5, columnspan=2)

        def get_out_entr():
            out_text_mg = ""
            sequence4 = str(seq_protentry.get(1.0, "end-1c"))
            sequence5 = sequence4.replace("\n", "")
            conc_uM = float(uM_conc_entr.get())
            if sequence5 != "":
                sequence6 = "".join([a for a in list(sequence5) if a in amino_acid_weights.keys()])
                mwt_calculated = 18.02 + sum([amino_acid_weights[a] for a in sequence6])
                res = conc_uM*float(mwt_calculated)/1000000
                out_text_mg += f"Concentration is {round(res, 2)} mg/ml, \nMolecular weight is {round(mwt_calculated, 2)} Da"
            elif sequence5 == "":
                mwt_provided = float(mwt_box_entr.get())
                res1 = conc_uM*float(mwt_provided)/1000000
                out_text_mg += f"Concentration is {round(res1, 2)} uM"
            seq_protout.delete('1.0', END)
            seq_protout.insert(END, out_text_mg)

        def clearall_3():
            seq_protout.delete('1.0', END)
            seq_protentry.delete('1.0', END)
            mwt_box_entr.delete('0', END)
            uM_conc_entr.delete('0', END)

if __name__ == "__main__":
    app = MultiPageApp()
    app.mainloop()
