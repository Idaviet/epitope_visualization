import pandas as pd
import re
from Bio.Data.IUPACData import protein_letters_1to3
from pymol import cmd

def create_epitope_selections_from_csv(epitope_csv, sequence, epitope_col='Epitope - Name', base_name='epitope'):
    """
    Creates a separate PyMOL selection for each row in an epitope CSV.

    For every row in the ``epitope_col`` column, extracts ``(AA, position)`` tuples
    via regex, converts each residue to its three-letter code, and builds a PyMOL
    selection spanning the range from the lowest to highest position. Each selection
    is displayed as a cyan cartoon.

    Args:
        epitope_csv (str): Path to the CSV file with an epitope column.
        sequence (str): Full antigen sequence (1-indexed). Note: not currently used
            for validation in this function.
        epitope_col (str): Column name containing residue strings (e.g.,
            ``"G446, Y449"``). Defaults to ``'Epitope - Name'``.
        base_name (str): Base name prefix for PyMOL selections. Defaults to
            ``'epitope'``. Note: the current implementation names selections as
            ``"{start}{aas}{end}"`` rather than using this parameter.

    Returns:
        None: Side effect only — creates and styles PyMOL selections.

    Example:
        Run from within the PyMOL command line after loading a structure::

            run pymol_functions.py
            create_epitope_selections_from_csv epitope_csv=/path/to/epitopes.csv, sequence="MYSFVSE..."
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
    """
    Creates a single PyMOL selection from all validated epitope residues across rows.

    Extracts all ``(AA, position)`` tuples from the CSV, validates each against the
    antigen sequence (1-based index, AA identity check), deduplicates, and builds one
    combined PyMOL selection. Residues that fall outside the sequence length are
    logged and skipped. Optionally scopes the selection to a specific chain.

    Args:
        epitope_csv (str): Path to the CSV file with an epitope column.
        sequence (str): Full antigen sequence used to validate each residue by position
            and identity.
        epitope_col (str): Column name containing residue strings (e.g.,
            ``"G446, Y449"``). Defaults to ``'Epitope - Name'``.
        selection_name (str): Name for the combined PyMOL selection. Defaults to
            ``'epitope_sel'``.
        chain (str or None): Optional chain ID (e.g., ``'A'``). If provided, clauses
            are scoped as ``"(chain {chain} and ...)"``. Defaults to ``None``.

    Returns:
        None: Side effect only — creates the selection and prints a summary.

    Raises:
        Prints a warning and returns early if no valid residues are found.

    Example:
        Run from within the PyMOL command line after loading a structure::

            run pymol_functions.py
            highlight_all_epitope_positions epitope_csv=/path/to/epitopes.csv, sequence="MYSFVSE...", chain="A"
    """
    
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
