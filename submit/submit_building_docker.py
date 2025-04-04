import argparse
import os
from tqdm import tqdm


INPUT_SMI_NAME = "input.smi"


SGE_TEMPLATE = """#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o {log_folder}
#$ -t 1-{count}
#$ -l h_rt={h_rt}
#$ -l mem_free=2.5G

export INDIR="{input_folder}"
{command}
"""

SLURM_TEMPLATE = """#!/bin/bash
#SBATCH --output={log_folder}/slurm-%A_%a.out
#SBATCH --array=1-{count}
#SBATCH --time={h_rt}
#SBATCH --mem=2500M

TMPDIR=$(mktemp -d /scratch/${{USER}}/job_${{SLURM_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}_XXXXXX)
trap "rm -rf $TMPDIR" EXIT
export INDIR="{input_folder}"
newgrp docker << EOF
{command}
EOF
"""

APPTAINER_COMMAND = "apptainer exec --cleanenv --no-mount tmp --bind ${{INDIR}}:/data --bind ${{TMPDIR}}:/tmp {container_path_or_name} bash /dock/ligand/submit/build-docker.sh"

DOCKER_COMMAND = "docker run --rm -u $(id -u):$(id -g) -v ${{INDIR}}:/data -v ${{TMPDIR}}:/tmp {container_path_or_name} bash /dock/ligand/submit/build-docker.sh"


def make_building_array_job(input_file, output_folder, bundle_size, minutes_per_mol, building_config_file,
                            array_job_name, skip_name_check, scheduler, container_software, container_path_or_name):
    all_ids = set()
    count = 1
    buffer = []
    os.makedirs(output_folder, exist_ok=True)
    with open(input_file) as f:
        for i,line in tqdm(enumerate(f), desc="Mols processed"):
            ll = line.split()
            if len(ll) != 2:
                raise ValueError(f"Input file {input_file} should have two columns: smiles and name (unique) in line {i+1}")
            smiles, name = ll
            if not skip_name_check and name in all_ids:
                raise ValueError(f"Name {name} is not unique.")
            if len(name) > 16:
                raise ValueError(f"Name {name} is too long. Max length is 16 characters.")
            if '.' in name:
                raise ValueError(f"Name {name} contains a period, which is not allowed.")
            if not skip_name_check:
                all_ids.add(name)
            buffer.append((smiles, name))
            if len(buffer) == bundle_size:
                output_one_list(buffer, count, output_folder)
                count += 1
                buffer = []
    if buffer:
        output_one_list(buffer, count, output_folder)
    else:
        count -= 1
    if count > 100000:
        raise ValueError(f"Too many molecules to build (array has {count} jobs, max is 100k). " +
                         "Increase bundle size or split the input file.")
    if count < 1000:
        print(f"WARNING: only {count} jobs in array. If you want your results fast, consider decreasing bundle size " +
              "to get more parallelization.")

    write_job_array_script(output_folder, count, bundle_size, minutes_per_mol, building_config_file, array_job_name, scheduler, container_software, container_path_or_name)

def write_job_array_script(output_folder, count, bundle_size, minutes_per_mol, building_config_file, 
                           array_job_name, scheduler, container_software, container_path_or_name):
    
    if building_config_file is not None:
        raise ValueError("Custom building config files are not yet supported")
    
    log_folder = os.path.join(output_folder, "logs")
    os.makedirs(log_folder, exist_ok=True)

    if container_software == "docker":
        command = DOCKER_COMMAND.format(
            container_path_or_name=container_path_or_name
        )
    elif container_software == "apptainer":
        command = APPTAINER_COMMAND.format(
            container_path_or_name=container_path_or_name
        )

    if scheduler == "sge":
        script = SGE_TEMPLATE.format(
            log_folder=log_folder,
            count=count,
            h_rt=minutes_to_h_rt(minutes_per_mol * bundle_size),
            input_folder=os.path.join(os.getcwd(),output_folder, "${SGE_TASK_ID}"),
            command=command,
        )
    elif scheduler == "slurm":
        script = SLURM_TEMPLATE.format(
            log_folder=log_folder,
            count=count,
            h_rt=minutes_to_h_rt(minutes_per_mol * bundle_size),
            input_folder=os.path.join(os.getcwd(),output_folder, "${SLURM_ARRAY_TASK_ID}"),
            command=command,
        )
    else:
        raise ValueError("Only SGE and Slurm are supported")

    with open(array_job_name, "w") as f:
        f.write(script)
    print(f"Array job script written to {array_job_name}")


def minutes_to_h_rt(minutes):
    if minutes > 60*24*14:
        raise ValueError("Requested time is too long (max is 14 days). Reduce bundle size or minutes_per_mol.")
    hours = minutes // 60
    remaining_minutes = int(minutes % 60)
    return f"{hours}:{remaining_minutes:02d}:00"


def output_one_list(buffer, count, output_folder):
    subfolder = os.path.join(output_folder, f"{count}")
    os.makedirs(subfolder, exist_ok=True)
    with open(os.path.join(subfolder, INPUT_SMI_NAME), "w") as f:
        for smiles, name in buffer:
            f.write(f"{smiles} {name}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Arguments for submitting a building array job from a single .smi file.")

    parser.add_argument("input_file", type=str, help="The input .smi file to build.")
    parser.add_argument("--scheduler", type=str, default='slurm', help="Which scheduler to use")
    parser.add_argument("--container_software", type=str, default="docker", help="Are we using docker or apptainer")
    parser.add_argument("--container_path_or_name", type=str, default="building_pipeline", help="Path to image for apptainer or name for docker")
    parser.add_argument("--output_folder", type=str, default='building_output',
                        help="The output folder to store the building results.")
    parser.add_argument("--bundle_size", type=int, default=1000,
                        help="The number of molecules to build per .db2.tgz bundle")
    parser.add_argument("--minutes_per_mol", type=float, default=3,
                        help="The time requested per molecule in minutes. Note that some molecules will have several " +
                        "protomers, so this should be well above the ~30 seconds per protomer that is typical.")
    parser.add_argument("--array_job_name", type=str, default="building_array_job.sh",
                        help="The name of the array job.")
    parser.add_argument("--building_config_file", type=str, help="Optional config file for building, to override " +
                        "default parameters.")
    parser.add_argument("--skip_name_check", action="store_true", help="If you know your molecule names are unique, skip the checks. " +
                        "This is useful for building lots of molecules where the set of names can take lots of memory")

    args = parser.parse_args()
    make_building_array_job(args.input_file, args.output_folder, args.bundle_size, args.minutes_per_mol,
                            args.building_config_file, args.array_job_name, args.skip_name_check, args.scheduler, 
                            args.container_software, args.container_path_or_name)


if __name__ == "__main__":
    main()
