#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import pandas as pd


TQ_VERBOSE = os.getenv("TQ_VERBOSE", "F").lower().startswith("t")
TQ_CACHE = os.getenv("TQ_CACHE", "F").lower().startswith("t")

class BLASTrunner: 

    def __init__(self, base_dir="/workspace/Test_Python_Jupyter/mi_entorno/docs", verbose=TQ_VERBOSE, tq_cache=TQ_CACHE):
        self.verbose = verbose
        self.tq_cache = tq_cache
        self.base_dir = base_dir
        self.cache_dir = os.path.join(self.base_dir, "blast_cache")
        self.output_dir_db = os.path.join(self.base_dir, "data_base")
        self.output_dir = os.path.join(self.base_dir, "blast_output")
        self.processed_sequences_file = os.path.join(self.cache_dir, "processed_sequences.fasta")

        os.makedirs(self.cache_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.output_dir_db, exist_ok=True)
    
    def _generate_cache_key(self, sequences):
        key = ''.join(sequences.keys())
        return hashlib.md5(key.encode()).hexdigest()

    def read_input(self, input_file):

        if self.verbose:
            print(f"Reading input file: {input_file}")

        try:
            with open(input_file, "r") as file:  
                lines = file.readlines()  

            sequences = {}
            sequence_name = None 

    
            for line in lines:
                line = line.strip()  
                if line.startswith(">"):  
                    sequence_name = line[1:] 
                    sequences[sequence_name] = ""  
                else:  
                    sequences[sequence_name] += line  

            return sequences  

        except FileNotFoundError:
            print(f"Error: The file {input_file} was not found.")  
            return None
        except Exception as e:
            print(f"Error: {e}")  
            return None

    def filter_new_and_removed_sequences(self, input_sequences):
        existing_sequences = {}
        if os.path.exists(self.processed_sequences_file):
            existing_sequences = self.read_input(self.processed_sequences_file)
        
        new_sequences = {k: v for k, v in input_sequences.items() if k not in existing_sequences}
        removed_sequences = {k: v for k, v in existing_sequences.items() if k not in input_sequences}
        
        return new_sequences, removed_sequences

    def update_cached_sequences(self, sequences):
        with open(self.processed_sequences_file, "w") as f:
            for name, seq in sequences.items():
                f.write(f">{name}\n{seq}\n")

    def make_blast_db(self, input_file, dbtype='prot', output_db='mi_base_de_datos'):
        input_sequences = self.read_input(input_file)
        if input_sequences is None:
            return

        new_sequences, removed_sequences = self.filter_new_and_removed_sequences(input_sequences)
        
        if not new_sequences and not removed_sequences:
            if self.verbose:
                print("No changes in sequences to update the BLAST database.")
            return

        if self.verbose:
            print(f"Updating BLAST database with {len(new_sequences)} new sequences and removing {len(removed_sequences)} sequences.")

        # Write all current sequences to a new file
        current_sequences_file = os.path.join(self.cache_dir, "current_sequences.fasta")
        with open(current_sequences_file, "w") as f:
            for name, seq in input_sequences.items():
                f.write(f">{name}\n{seq}\n")

        command = [
            "makeblastdb",
            "-in", current_sequences_file,
            "-dbtype", dbtype,
            "-out", os.path.join(self.output_dir_db, output_db)
        ]
        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode == 0:
            self.update_cached_sequences(input_sequences)
            if self.verbose:
                print("BLAST database updated successfully.")
        else:
            print("Error updating BLAST database.", result.stderr)


    def run_blastp(self, input_file, database_prefix, output_file):

        if self.verbose:
            print(f"Running BLASTP with input file: {input_file} and database: {database_prefix}")

        command = [
            "blastp",
            "-query", input_file,
            "-db", database_prefix,
            "-out", output_file,
            "-outfmt", "6"  # Formato tabular
        ]
        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode == 0:
            if self.verbose:
                print(f"BLASTP completed successfully. Results saved to {output_file}.")
        else:
            print("Error running BLASTP.", result.stderr)

if __name__ == "__main__":
    BLAST = BLASTrunner(verbose=TQ_VERBOSE, tq_cache=TQ_CACHE)
    
    input_file = "/workspace/Test_Python_Jupyter/mi_entorno/docs/GH51_short.txt"
    database_prefix = "/workspace/Test_Python_Jupyter/mi_entorno/docs/data_base/mi_base_de_datos"
    output_file = "/workspace/Test_Python_Jupyter/mi_entorno/docs/blast_output/blast_results.txt"
    
    # Leer el archivo de entrada
    sequences = BLAST.read_input(input_file)
    print(len(sequences))

    # Crear la base de datos BLAST
    BLAST.make_blast_db(input_file)

    # Ejecutar BLASTP
    BLAST.run_blastp(input_file, database_prefix, output_file)