#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the GSE25066 analysis pipeline.")
    parser.add_argument(
        "--config",
        default="config/analysis_config.yaml",
        help="Path to the pipeline config file.",
    )
    parser.add_argument(
        "--skip-modeling",
        action="store_true",
        help="Run preprocessing/EDA/BI/summary only and skip Python modeling.",
    )
    return parser.parse_args()


def run_command(cmd: list[str], cwd: Path, step_name: str) -> None:
    print(f"[pipeline] start: {step_name}")
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        raise SystemExit(f"[pipeline] failed: {step_name} (exit code {result.returncode})")
    print(f"[pipeline] done: {step_name}")


def ensure_conda_env(root_dir: Path, env_name: str) -> None:
    conda_bin = shutil.which("conda")
    if conda_bin is None:
        raise SystemExit("conda was not found in PATH.")

    result = subprocess.run(
        [conda_bin, "env", "list"],
        cwd=root_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise SystemExit("failed to inspect conda environments.")

    env_names = {line.split()[0] for line in result.stdout.splitlines() if line and not line.startswith("#")}
    if env_name not in env_names:
        raise SystemExit(
            f"conda environment '{env_name}' was not found. "
            f"Create it with: conda env create -n {env_name} -f config/bc_environment.yml"
        )


def main() -> None:
    args = parse_args()
    root_dir = Path(__file__).resolve().parent
    config_path = (root_dir / args.config).resolve() if not Path(args.config).is_absolute() else Path(args.config)

    if not config_path.exists():
        raise SystemExit(f"config file not found: {config_path}")

    r_steps = [
        ("prepare_data", ["Rscript", str(root_dir / "scripts/R/01_prepare_data.R"), str(config_path)]),
        ("eda", ["Rscript", str(root_dir / "scripts/R/02_eda.R"), str(config_path)]),
        ("bi_analysis", ["Rscript", str(root_dir / "scripts/R/03_bi_analysis.R"), str(config_path)]),
    ]

    for step_name, cmd in r_steps:
        run_command(cmd, cwd=root_dir, step_name=step_name)

    if not args.skip_modeling:
        ensure_conda_env(root_dir, "bc")
        run_command(
            [
                "conda",
                "run",
                "-n",
                "bc",
                "python",
                str(root_dir / "scripts/python/train_models.py"),
                "--config",
                str(config_path),
            ],
            cwd=root_dir,
            step_name="modeling",
        )

    run_command(
        ["Rscript", str(root_dir / "scripts/R/04_compile_summary.R"), str(config_path)],
        cwd=root_dir,
        step_name="compile_summary",
    )

    print("[pipeline] completed")


if __name__ == "__main__":
    main()
