"""Command line interface:  python -m thermocalculus <command> ..."""
import argparse
import sys

from . import __version__, scenario, symbolic
from .eos import load_fluids
from .errors import ThermoError


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python -m thermocalculus",
        description="Computational engine for the calculus of thermodynamics")
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run one or more scenario files")
    r.add_argument("files", nargs="+")
    r.add_argument("-o", "--outdir", default="outputs")
    d = sub.add_parser("derive", help="reduce a derivative, e.g. \"(dT/dP)_H\"")
    d.add_argument("derivative")
    sub.add_parser("maxwell", help="verify the four Maxwell relations")
    sub.add_parser("fluids", help="list the fluid database")
    a = sub.add_parser("ask", help="interactive session: the engine asks for missing data")
    a.add_argument("--script", help="file of pre-written answers instead of the keyboard")
    a.add_argument("--fluid", help="fluid for an EoS model, e.g. CO2")
    a.add_argument("--eos", default="PR", help="IG | vdW | RK | SRK | PR (default PR)")
    a.add_argument("--transcript", help="also write the session to this text file")
    args = ap.parse_args(argv)

    if args.cmd == "run":
        worst = 0
        for path in args.files:
            res = scenario.run(path)
            out_csv, out_log = scenario.write(res, args.outdir, path)
            print(f"[{res.status}] {res.name} -> {out_csv}")
            for line in res.log:
                if line.split(":")[0] in ("ERROR", "MISSING DATA", "INCONSISTENT",
                                          "WARNING", "SINGULAR", "NOTE") or line.startswith("ERROR"):
                    print("    " + line)
            if res.status.startswith("error"):
                worst = 1
        return worst
    if args.cmd == "derive":
        try:
            expr = symbolic.reduce_derivative(args.derivative)
        except ThermoError as exc:
            print("error:", exc)
            return 1
        print(f"{args.derivative} = {expr}")
        print("needs:", ", ".join(sorted(expr.symbols())) or "nothing")
        print("LaTeX:", expr.latex())
        return 0
    if args.cmd == "maxwell":
        for text, origin, value, ok in symbolic.check_maxwell():
            print(f"{'ok  ' if ok else 'FAIL'} {text} = {value}    (from {origin})")
        return 0
    if args.cmd == "ask":
        from . import interactive
        lines = []

        def out(text):
            print(text)
            lines.append(text)
        interactive.run_session(args.script, args.fluid, args.eos, out)
        if args.transcript:
            with open(args.transcript, "w") as fh:
                fh.write("\n".join(lines) + "\n")
        return 0
    if args.cmd == "fluids":
        for f in load_fluids().values():
            cp = (f"Cp0 {f.cp0.tmin:g}-{f.cp0.tmax:g} K" if f.cp0 else "no Cp0 model")
            print(f"{f.name:16s} Tc = {f.Tc:9.4f} K  Pc = {f.Pc / 1e6:7.4f} MPa  "
                  f"omega = {f.omega:7.5f}  M = {f.M * 1e3:8.4f} g/mol  {cp}")
        return 0


if __name__ == "__main__":
    sys.exit(main())
