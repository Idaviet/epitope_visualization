import pandas as pd
import re
from Bio.Data.IUPACData import protein_letters_1to3
from pymol import cmd

def create_epitope_selections_from_csv(epitope_csv, sequence, epitope_col='Epitope - Name', base_name='epitope'):
    """
    For each row of epitopes in the CSV, create a separate PyMOL selection
    with validated residues based on sequence identity and position.
    
    :param epitope_csv: Path to the CSV with an 'epitope' column
    :param sequence: Full antigen sequence (1-indexed)
    :param chain: Chain ID to use in PyMOL selections
    :param epitope_col: Column in CSV containing residue strings (e.g. "G446, Y449")
    :param base_name: Base name for PyMOL selections (will be suffixed with row number)
    """
    df = pd.read_csv(epitope_csv)
    
    for i, row in enumerate(df[epitope_col].dropna()):
        aas = ""
        poss = []
        residues = re.findall(r"([A-Z])(\d+)", row)
        clauses = []

        for aa, pos in residues:
            aas += aa
            pos = int(pos)
            poss.append(pos)
            aa3 = protein_letters_1to3.get(aa.upper(), "UNK")
            clause = f"(resi {pos} and resn {aa3})"
            clauses.append(clause)
        
        if clauses:
            start = min(poss)
            end = max(poss)

            selection_str = " or ".join(clauses)
            selection_name = f"{start}{aas}{end}"
            cmd.select(selection_name, selection_str)
            cmd.show("cartoon", selection_name)
            cmd.color("cyan", selection_name)

def highlight_all_epitope_positions(epitope_csv, sequence, epitope_col="Epitope - Name", selection_name="epitope_sel", chain=None):
    df = pd.read_csv(epitope_csv)
    matched = set()

    for row in df[epitope_col].dropna():
        residues = re.findall(r"([A-Z])(\d+)", row)
        for aa, pos in residues:
            pos = int(pos)
            if 1 <= pos <= len(sequence):
                if sequence[pos - 1] == aa:
                    matched.add((aa, pos))
            else:
                print(f"[!] Position {pos} out of range")

    if not matched:
        print("[!] No valid epitope positions found.")
        return

    sel_clauses = []
    for aa, pos in sorted(matched, key=lambda x: x[1]):
        aa3 = protein_letters_1to3.get(aa.upper(), "UNK")
        clause = f"(resi {pos} and resn {aa3})"
        if chain:
            clause = f"(chain {chain} and {clause})"
        sel_clauses.append(clause)

    selection_str = " or ".join(sel_clauses)
    cmd.select(selection_name, selection_str)
    print(f"Selection '{selection_name}' created with {len(matched)} residues.")
