
from __future__ import annotations
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="microctpm command-line tools")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("rev-xy")
    p.add_argument("image"); p.add_argument("output_dir"); p.add_argument("--voxel-size", type=float, default=46); p.add_argument("--threshold", type=float, default=55)
    p = sub.add_parser("rev-xyz")
    p.add_argument("image"); p.add_argument("output_dir"); p.add_argument("--voxel-size", type=float, default=15); p.add_argument("--threshold", type=float, default=85)
    p = sub.add_parser("pore-metrics")
    p.add_argument("image"); p.add_argument("output_dir"); p.add_argument("--voxel-size", type=float, default=15); p.add_argument("--threshold", type=float, default=85)
    p = sub.add_parser("critical-pore")
    p.add_argument("image"); p.add_argument("output_dir"); p.add_argument("--voxel-size", type=float, default=15e-6); p.add_argument("--threshold", type=float, default=110); p.add_argument("--r-max", type=int, default=10); p.add_argument("--sigma", type=float, default=1)
    p = sub.add_parser("psd-wrc")
    p.add_argument("image"); p.add_argument("output_dir"); p.add_argument("--voxel-size", type=float, default=46); p.add_argument("--threshold", type=float, default=85)
    p = sub.add_parser("permeability")
    p.add_argument("image"); p.add_argument("output_dir"); p.add_argument("--voxel-size", type=float, default=15e-6); p.add_argument("--threshold", type=float, default=85)

    args = parser.parse_args()
    if args.command == "rev-xy":
        from .rev_xy import run_rev_xy; run_rev_xy(args.image, args.output_dir, args.voxel_size, args.threshold)
    elif args.command == "rev-xyz":
        from .rev_xyz import run_rev_xyz; run_rev_xyz(args.image, args.output_dir, args.voxel_size, args.threshold)
    elif args.command == "pore-metrics":
        from .pore_metrics import run_pore_metrics; run_pore_metrics(args.image, args.output_dir, args.voxel_size, args.threshold)
    elif args.command == "critical-pore":
        from .critical_pore import calculate_critical_pore_diameter; calculate_critical_pore_diameter(args.image, args.voxel_size, args.threshold, args.r_max, args.sigma, args.output_dir)
    elif args.command == "psd-wrc":
        from .pore_size_distribution import run_pore_size_distribution; run_pore_size_distribution(args.image, args.output_dir, args.voxel_size, args.threshold)
    elif args.command == "permeability":
        from .permeability import calculate_permeability; calculate_permeability(args.image, args.output_dir, args.voxel_size, args.threshold)

if __name__ == "__main__":
    main()
