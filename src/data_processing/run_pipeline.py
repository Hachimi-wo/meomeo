import os
import subprocess
import sys

def run_script(script_path):
    print(f"\n{'='*60}\nRunning: {script_path}\n{'='*60}")
    result = subprocess.run([sys.executable, script_path], cwd=os.path.dirname(script_path))
    if result.returncode != 0:
        print(f"Error executing {script_path}")
        sys.exit(1)

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    phase2_dir = os.path.join(base_dir, "01_build_dataset")
    
    # Check if 07 output already exists to populate cautious cache for 06b
    out_07 = os.path.abspath(os.path.join(phase2_dir, "../../../data/01_intermediate/phase2_assembly/07_pubchem_scored_negatives.csv"))
    cache_file = os.path.abspath(os.path.join(phase2_dir, "../../../data/01_intermediate/phase2_assembly/cache_pubchem_scores.csv"))
    
    if os.path.exists(out_07) and not os.path.exists(cache_file):
        print(f"Populating Cautious Cache from existing output: {out_07}")
        import shutil
        shutil.copy(out_07, cache_file)
    
    scripts = [
        "01_merge_adrecs_sider.py",
        "02_annotate_labels_from_pattern.py",
        "03_deduplicate_by_drugname.py",
        "04_augment_drugcentral_smiles.py",
        "05_integrate_ctd_markers.py",
        "06a_create_candidate_negatives.py",
        "06b_score_candidates_cautious_cache.py",  # Use the new cautious cache version
        "07_assemble_final_dataset.py"
    ]
    
    for script in scripts:
        script_path = os.path.join(phase2_dir, script)
        run_script(script_path)
        
    print("\n[SUCCESS] Pipeline 01_build_dataset completed successfully!")

if __name__ == "__main__":
    main()
