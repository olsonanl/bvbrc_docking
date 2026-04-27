#
# Support code for diffdock v1.1
#
# Here we don't need to separately create the encodings as it is done internally
# by diffdock. However the inference run must be performed from
# the diffdock build/install directory. When using the
# BV-BRC container the location of this directory may be found
# in the BVDOCK_DIFFDOCK_DIR environment variable.
#


import os
import re
from operator import itemgetter

import numpy as np
import pandas as pd
from tqdm import tqdm

from bvbrc_docking.utils import (
    cal_cnn_aff,
    clean_pdb,
    comb_pdb,
    run_and_save,
    run_list_and_save,
    sdf2pdb,
    validate_smiles,
)


class diff_dock(object):
    def __init__(
        self,
        receptor_pdb,
        drug_dbs,
        diffdock_dir,
        output_dir,
        top_n: int = 1,
        batch_size: int = -1,
        cont_run=False,
        **kwargs,
    ) -> None:
        self.receptor_pdb = os.path.abspath(receptor_pdb)
        self.label = os.path.basename(receptor_pdb).split(".")[0]
        self.drug_dbs = os.path.abspath(drug_dbs)
        self.diffdock_dir = os.path.abspath(diffdock_dir)
        self.output_dir = os.path.abspath(output_dir)
        if cont_run:
            os.makedirs(self.output_dir, exist_ok=True)
        else:
            os.makedirs(self.output_dir)

        self.run_dir = self.output_dir

        self.top_n = top_n
        self.batch_size = batch_size

    def prepare_inputs(self):
        inputs = []
        failed = 0
        with open(self.drug_dbs, "r") as fp:
            for line in fp:
                # add support for single col of smile strings HM
                ident, smiles_str = line.split()
                if not validate_smiles(smiles_str):
                    #
                    # See if fields were reversed
                    #
                    if validate_smiles(ident):
                        ident, smiles_str = smiles_str, ident
                    else:
                        failed += 1
                        print(f"Smiles string for compound f{ident} is not valid")
                        continue
                inputs.append([ident, smiles_str])
        if failed:
            print(
                f"Failure f{failed} parsing smiles strings from input file f{self.drug_dbs}"
            )
            os.exit(1)

        output_pdb = os.path.join(self.run_dir, f"{self.label}.pdb")
        self.pdb_file = clean_pdb(self.receptor_pdb, output_pdb)

        self.all_runs = f"{self.run_dir}/all.csv"
        with open(self.all_runs, "w") as out:
            out.write("protein_path,ligand_description,complex_name,protein_sequence\n")
            for ident, smiles_str in inputs:
                out.write(f"{self.pdb_file},{smiles_str},{ident}\n")
        return inputs

    def run(self):
        log_file = f"{self.run_dir}/diffdock_log"
        self.log_handle = open(log_file, "w")

        input_set = self.prepare_inputs()
        self.run_docking()
        self.post_process(input_set)

        self.log_handle.close()

    def run_docking(self):
        #
        # Run with -u to enable line-buffering so we can see the errors aligned with
        # the stdout for better debugging.
        #
        cmd_diffdock = [
            "python",
            "-u",
            "-m", "inference",
            "--protein_ligand_csv", self.all_runs,
            "--out_dir", self.run_dir,
            ]

        if self.batch_size > 0:
            cmd_diffdock.extend(["--batch_size", str(self.batch_size)])
            
        # the run failed with the original diffdock 1.0 parameters used:
        # f"--inference_steps 20 --samples_per_complex 40 --batch_size 6"

        proc = run_list_and_save(
            cmd_diffdock, cwd=self.diffdock_dir, output_file=self.log_handle
        )

    def post_process(self, input_set):
        #
        # Results are in directories named by the identifiers
        #

        for ident, smiles_str in tqdm(input_set):
            by_rank = []
            result_path = f"{self.run_dir}/{ident}"

            for file in os.listdir(result_path):
                m = re.match(r"rank(\d+)_confidence(-?\d\.\d+).sdf", file)
                if m:
                    rank, confidence = m.group(1, 2)

                    # skip high score and rank
                    rank = int(rank)
                    if self.top_n != 0 and rank > self.top_n:
                        continue
                    if float(confidence) > 100:
                        continue

                    # convert ligand sdf to complex with protein
                    sdf_file = os.path.join(result_path, file)
                    lig_pdb = sdf2pdb(sdf_file)
                    combined = comb_pdb(self.pdb_file, lig_pdb)
                    if combined is None:
                        continue
                    else:
                        by_rank.insert(
                            0,
                            [
                                ident,
                                os.path.join(result_path, file),
                                rank,
                                confidence,
                                combined,
                            ],
                        )

            by_rank.sort(key=itemgetter(2))

            # for idx, ent in enumerate(by_rank):
            #     print(f"idx={idx} ent={ent}")
            #     ident, file, rank, confidence = ent

            #     #
            #     # For the final output, we combine each of the
            #     # docked ligand with the original PDF for easy viewing
            #     #
            #     # We need to convert the sdf to pdb first.
            #     #
            #     lig_pdb = sdf2pdb(file)

            #     combined = comb_pdb(self.pdb_file, lig_pdb)
            #     if combined is None:
            #         by_rank.remove(ent)
            #         continue
            #     ent.append(combined)

            with open(f"{result_path}/result.csv", "w") as fp:
                print(
                    "\t".join(
                        [
                            "ident",
                            "rank",
                            "score",
                            "lig_sdf",
                            "comb_pdb",
                            "CNNscore",
                            "CNNaffinity",
                            "Vinardo",
                        ]
                    ),
                    file=fp,
                )
                for ident, path, rank, confidence, combined_path in by_rank:
                    mol = cal_cnn_aff(
                        self.pdb_file,
                        path,
                        f"{self.diffdock_dir}/gnina",
                        log_handle=self.log_handle,
                    )
                    if mol is None:
                        continue
                    print(
                        "\t".join(
                            [
                                ident,
                                str(rank),
                                str(confidence),
                                os.path.basename(path),
                                os.path.basename(combined_path),
                                str(mol.data["CNNscore"]),
                                str(mol.data["CNNaffinity"]),
                                str(mol.data["minimizedAffinity"]),
                            ]
                        ),
                        file=fp,
                    )

        #     for j in range(self.top_n):
        #         sdf = glob.glob(
        #             f"{glob.escape(result_path)}/rank{j+1}_confidence-*.sdf"
        #         )[0]
        #         score = float(os.path.basename(sdf).split("-")[1][:-4])
        #         local_dict = row.to_dict()
        #         local_dict["lig_sdf"] = sdf
        #         local_dict["score"] = score
        #         local_dict["comp_pdb"] = comb_pdb(
        #             local_dict["protein_path"], sdf2pdb(sdf)
        #         )
        #         output_df.append(local_dict)

        # output_df = pd.DataFrame(output_df)
        # output_df.to_csv(f"{self.run_dir}/result.csv")
        # return output_df
