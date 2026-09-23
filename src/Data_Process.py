import os
import subprocess
import sys
import time

def run_script(script_path, output_check=None):
    script_name = os.path.basename(script_path)
    if output_check and os.path.exists(output_check):
        print(f"[CACHE]   {script_name:<40} -> Skipped")
        return

    print(f"[RUN  ]   {script_name:<40} ", end="", flush=True)
    
    result = subprocess.run(
        [sys.executable, script_path], 
        cwd=os.path.dirname(script_path), 
        capture_output=True, 
        text=True
    )
    
    if result.returncode != 0:
        print("[FAIL ]")
        print(f"\n{'-'*60}\nERROR LOG FOR: {script_name}\n{'-'*60}")
        print(result.stdout)
        print(result.stderr)
        print('-'*60)
        sys.exit(1)
    else:
        print("[OK   ]")

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    
    # CREATE REQUIRED DIRECTORIES
    os.makedirs(os.path.join(base_dir, "../data/02_curated_raw"), exist_ok=True)
    os.makedirs(os.path.join(base_dir, "../data/03_final_dataset"), exist_ok=True)
    os.makedirs(os.path.join(base_dir, "../data/03_precomputed_features"), exist_ok=True)
    
    # Intermediate subdirs
    for phase in ["phase1_curation", "phase2_assembly"]:
        os.makedirs(os.path.join(base_dir, f"../data/01_intermediate/{phase}"), exist_ok=True)
    for sub in ["adrecs", "sider", "drugcentral", "ctdbase"]:
        os.makedirs(os.path.join(base_dir, f"../data/01_intermediate/phase1_curation/{sub}"), exist_ok=True)

    # PHASE 1: RAW TO CURATED (00_raw_to_curated)
    p1 = os.path.join(base_dir, "data_processing", "00_raw_to_curated")
    
    # Adrecs
    d_adrecs = os.path.join(p1, "01_adrecs")
    o_adrecs = "../../../../data/01_intermediate/phase1_curation/adrecs"
    run_script(os.path.join(d_adrecs, "01_adrecs_merge.py"), os.path.join(d_adrecs, f"{o_adrecs}/merged_database_v2.csv"))
    run_script(os.path.join(d_adrecs, "02_adrecs_fetch_pubchem.py"), os.path.join(d_adrecs, f"{o_adrecs}/adrecs_pubchem_smiles.csv"))
    run_script(os.path.join(d_adrecs, "03_adrecs_fetch_chembl.py"), os.path.join(d_adrecs, f"{o_adrecs}/adrecs_chembl_smiles.csv"))
    run_script(os.path.join(d_adrecs, "04_adrecs_run_mehc.py"), os.path.join(d_adrecs, f"{o_adrecs}/mehc/refinement/post_refined_smiles.csv"))
    run_script(os.path.join(d_adrecs, "05_adrecs_merge_metadata.py"), os.path.join(d_adrecs, "../../../../data/02_curated_raw/adrecs_pneumotox_curated.csv"))

    # Sider
    d_sider = os.path.join(p1, "02_sider")
    o_sider = "../../../../data/01_intermediate/phase1_curation/sider"
    run_script(os.path.join(d_sider, "01_sider_fetch_pubchem.py"), os.path.join(d_sider, f"{o_sider}/sider_filtered_nomehc.csv"))

    run_script(os.path.join(d_sider, "03_sider_prepare_mehc.py"), os.path.join(d_sider, f"{o_sider}/respiratory_labeled_for_curation.csv"))
    run_script(os.path.join(d_sider, "04_sider_run_mehc.py"), os.path.join(d_sider, f"{o_sider}/mehc/refinement/post_refined_smiles.csv"))
    run_script(os.path.join(d_sider, "05_sider_merge_metadata.py"), os.path.join(d_sider, "../../../../data/02_curated_raw/sider_respiratory_adrs.csv"))

    # DrugCentral
    d_dc = os.path.join(p1, "03_drugcentral")
    o_dc = "../../../../data/01_intermediate/phase1_curation/drugcentral"
    run_script(os.path.join(d_dc, "01_drugcentral_filter.py"), os.path.join(d_dc, f"{o_dc}/drugcentralonly.csv"))
    run_script(os.path.join(d_dc, "02_drugcentral_run_mehc.py"), os.path.join(d_dc, "../../../../data/02_curated_raw/drugcentral_smiles.csv"))

    # CTDBase
    d_ctd = os.path.join(p1, "04_ctdbase")
    o_ctd = "../../../../data/01_intermediate/phase1_curation/ctdbase"
    run_script(os.path.join(d_ctd, "01_ctd_fetch_api.py"), os.path.join(d_ctd, f"{o_ctd}/ctdbase_smiles.csv"))
    run_script(os.path.join(d_ctd, "02_ctd_run_mehc.py"), os.path.join(d_ctd, f"{o_ctd}/mehc/refinement/post_refined_smiles.csv"))
    run_script(os.path.join(d_ctd, "03_ctd_merge_metadata.py"), os.path.join(d_ctd, "../../../../data/02_curated_raw/ctd_mechanism_markers.csv"))

    # PHASE 2: BUILD DATASET (01_build_dataset)
    p2 = os.path.join(base_dir, "data_processing", "01_build_dataset")
    
    # Populate cautious cache for 06b if needed
    out_07 = os.path.join(p2, "../../../data/01_intermediate/phase2_assembly/07_pubchem_scored_negatives.csv")
    cache_file = os.path.join(p2, "../../../data/01_intermediate/phase2_assembly/cache_pubchem_scores.csv")
    if os.path.exists(out_07) and not os.path.exists(cache_file):
        print(f"Cautious Cache from existing output: {out_07}")
        import shutil
        shutil.copy(out_07, cache_file)
        
    run_script(os.path.join(p2, "01_merge_adrecs_sider.py"))
    run_script(os.path.join(p2, "02_annotate_labels_from_pattern.py"))
    run_script(os.path.join(p2, "03_deduplicate_by_drugname.py"))
    run_script(os.path.join(p2, "04_augment_drugcentral_smiles.py"))
    run_script(os.path.join(p2, "05_integrate_ctd_markers.py"))
    run_script(os.path.join(p2, "06a_create_candidate_negatives.py"))
    run_script(os.path.join(p2, "06b_score_candidates_cautious_cache.py"))
    run_script(os.path.join(p2, "07_assemble_final_dataset.py"))
    
    print("\n done ")

if __name__ == "__main__":
    main()
