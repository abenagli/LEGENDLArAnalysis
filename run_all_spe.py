#!/usr/bin/env python3

import json
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass

SCRIPT_DIR = Path(__file__).parent.resolve()
CONFIG_FILE = SCRIPT_DIR / "config.json"
MACRO_FILE = SCRIPT_DIR / "adc_read_all_data.C"


@dataclass
class Job:
    dataset: str
    run: int
    charge_time: int
    data_dir: Path
    output_dir: Path


def load_config():
    """Load configuration from config.json"""
    with open(CONFIG_FILE, "r") as f:
        return json.load(f)


def generate_jobs(config):
    """Generate all processing jobs from config"""
    jobs = []
    base_path = Path(config["base_path"])
    
    for dataset in config["datasets"]:
        dataset_name = dataset["name"]
        runs = dataset["runs"]
        charge_times = dataset["charge_times"]
        
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
        f'{MACRO_FILE}("{job.data_dir}/",{job.charge_time})'
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
    
    config = load_config()
    jobs = generate_jobs(config)
    n_workers = config.get("n_workers", 8)
    
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

