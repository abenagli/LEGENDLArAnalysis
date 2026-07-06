#!/usr/bin/env python3

import shutil
import subprocess
import argparse
import importlib.util
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass

SCRIPT_DIR = Path(__file__).parent.resolve()
MACRO_FILE = SCRIPT_DIR / "adc_read_all_data.C"


@dataclass
class Job:
    dataset: str
    run: int
    charge_time: int
    start_charge: float
    data_dir: Path
    output_dir: Path


def load_config(config_file):
    """Load configuration from Python config file"""
    config_path = Path(config_file)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    # Load Python module dynamically
    spec = importlib.util.spec_from_file_location("config_module", config_path)
    config_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(config_module)
    
    return {
        "base_path": config_module.BASE_PATH,
        "n_workers": config_module.N_WORKERS,
        "macro_file": config_module.MACRO_FILE,
        "datasets": config_module.CONFIG["datasets"]
    }


def generate_jobs(config):
    """Generate all processing jobs from config"""
    jobs = []
    base_path = Path(config["base_path"])
    
    for dataset in config["datasets"]:
        dataset_name = dataset["name"]
        runs = dataset["runs"]
        charge_times = dataset["charge_times"]
        start_charge = dataset.get("startCharge", 4000)  # default to 4000 if not specified
        
        for run in runs:
            run_dir = base_path / dataset_name / f"run{run}"
            
            if not run_dir.exists():
                print(f"WARNING: {run_dir} does not exist, skipping")
                continue
            
            for charge_time in charge_times:
                job = Job(
                    dataset=dataset_name,
                    run=run,
                    charge_time=charge_time,
                    start_charge=start_charge,
                    data_dir=run_dir,
                    output_dir=run_dir
                )
                jobs.append(job)
    
    return jobs


def run_one(job):
    """Process a single job"""
    workdir = SCRIPT_DIR / f"tmp_job_{job.dataset}_run{job.run}_{job.charge_time}ns"
    
    if workdir.exists():
        shutil.rmtree(workdir)
    
    workdir.mkdir()
    
    cmd = [
        "root",
        "-l",
        "-b",
        "-q",
        f'{MACRO_FILE}("{job.data_dir}/",{job.charge_time},{job.start_charge})'
    ]
    
    label = f"[{job.dataset}/run{job.run}] {job.charge_time} ns"
    print(f"{label} starting", flush=True)
    
    result = subprocess.run(
        cmd,
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    
    produced_file = workdir / "analyzed.root"
    final_name = job.output_dir / f"analyzed_{job.charge_time}ns.root"
    
    if result.returncode != 0:
        shutil.rmtree(workdir, ignore_errors=True)
        return job, False, f"ROOT execution failed:\n{result.stdout}"
    
    if not produced_file.exists():
        shutil.rmtree(workdir, ignore_errors=True)
        return job, False, f"analyzed.root not found:\n{result.stdout}"
    
    if final_name.exists():
        final_name.unlink()
    
    shutil.move(str(produced_file), str(final_name))
    shutil.rmtree(workdir, ignore_errors=True)
    
    return job, True, str(final_name)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description="Process ADC data with configurable datasets and parameters",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                        # Use default config.py in current directory
  %(prog)s config.py              # Specify config file
  %(prog)s /path/to/config.py     # Use absolute path
        """
    )
    
    parser.add_argument(
        "config",
        nargs="?",
        default=str(SCRIPT_DIR / "config.py"),
        help="Path to config.py file (default: ./config.py)"
    )
    
    args = parser.parse_args()
    
    try:
        config = load_config(args.config)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        exit(1)
    except Exception as e:
        print(f"ERROR: Failed to load config: {e}")
        exit(1)
    
    jobs = generate_jobs(config)
    n_workers = config.get("n_workers", 8)
    
    if not jobs:
        print("No jobs to process (check your config file)")
        exit(0)
    
    print(f"Found {len(jobs)} jobs to process\n")
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        
        futures = {
            executor.submit(run_one, job): job
            for job in jobs
        }
        
        completed = 0
        failed = 0
        
        for future in as_completed(futures):
            job, ok, message = future.result()
            
            label = f"[{job.dataset}/run{job.run}] {job.charge_time} ns"
            
            if ok:
                print(f"{label} ✓ done")
                completed += 1
            else:
                print(f"{label} ✗ FAILED")
                print(f"  {message}\n")
                failed += 1
    
    print(f"\n=== Summary ===")
    print(f"Completed: {completed}/{len(jobs)}")
    print(f"Failed: {failed}/{len(jobs)}")

