"""Interactive session: the engine asks the user for the data it is missing.

Topic 2 requires software that "identifies missing data, prompts the user if
additional information is required, and performs the computation".  The
scenario files do the first and last part in batch mode; this module adds the
prompting loop.

Session outline
---------------
1. The user states the targets: a derivative such as ``(dT/dP)_H`` or a
   quantity such as ``Cv``, ``w`` or ``muJT``.
2. The user enters what is known, one ``name = value unit`` per line.
3. The engine forward-chains the exact relations of ``knowledge.py`` (and,
   when a fluid and an equation of state are chosen with ``model``, evaluates
   the EoS as soon as T and P are known).
4. For every target that cannot yet be computed, the engine lists the smallest
   sets of additional quantities that would make it computable and prompts for
   one of them.  The answer is added and the analysis is repeated.
5. When nothing more can be learned, the results are printed with provenance.

``input_fn`` and ``output_fn`` are injectable, so a session can be driven from
a script file (``python -m thermocalculus ask --script file``) and recorded.
"""
from . import knowledge, props
from .eos import EOS_NAMES, EoS, get_fluid
from .errors import ThermoError
from .symbolic import parse_derivative
from .units import SI_UNIT, parse_quantity

HELP = """\
Commands (type one at the prompt):
  <name> = <value> <unit>   state a known quantity, e.g.  T = 350 K  or  kappaT = 0.15 1/MPa
  target <item>             add a target: a quantity (Cv, w, muJT, ...) or a derivative (dT/dP)_H
  model <fluid> <eos>       evaluate an equation of state when T and P are known, e.g. model CO2 PR
  show                      list what is known and what has been inferred
  go                        analyse now (also done automatically after every input)
  help                      this text
  quit                      end the session
Quantities: T P V Z M rho S Cp Cv gamma alpha kappaT kappaS gammaV w muJT"""


class Session:
    def __init__(self, input_fn=input, output_fn=print):
        self.ask_user = input_fn
        self.say = output_fn
        self.given = {}          # name -> SI value, stated by the user
        self.targets = []        # quantity names or derivative strings
        self.model = None        # (fluid, eos) or None
        self.model_name = ""
        self.skipped = set()     # targets the user declined to pursue

    # ---- input handling -------------------------------------------------
    def read_line(self, prompt):
        try:
            return self.ask_user(prompt)
        except EOFError:
            return None

    def add_target(self, item):
        item = item.strip()
        if not item:
            return
        if item in knowledge.QUANTITIES:
            ok = True
        else:
            try:
                parse_derivative(item)
                ok = True
            except ThermoError as exc:
                self.say(f"  cannot use '{item}' as a target: {exc}")
                return
        if ok and item not in self.targets:
            self.targets.append(item)
            self.say(f"  target added: {item}")

    def add_value(self, line):
        name, _, text = line.partition("=")
        name, text = name.strip(), text.strip()
        if name not in knowledge.QUANTITIES:
            self.say(f"  '{name}' is not a quantity the engine knows. "
                     f"Known: {' '.join(knowledge.QUANTITIES)}")
            return False
        try:
            self.given[name] = parse_quantity(text, name)
        except ThermoError as exc:
            self.say(f"  rejected: {exc}")
            return False
        self.say(f"  {name} = {text}   (stored as {self.given[name]:.6g} {SI_UNIT.get(name, '')})")
        return True

    def set_model(self, text):
        parts = text.split()
        if len(parts) != 2:
            self.say("  usage: model <fluid> <eos>   (eos: " + ", ".join(EOS_NAMES) + ")")
            return
        try:
            fluid = get_fluid(parts[0])
            eos = EoS(parts[1], fluid.Tc, fluid.Pc, fluid.omega)
        except ThermoError as exc:
            self.say(f"  rejected: {exc}")
            return
        self.model, self.model_name = (fluid, eos), f"{fluid.name} / {parts[1]}"
        self.say(f"  model set: {self.model_name}. T and P are now enough to evaluate the state.")

    # ---- analysis -------------------------------------------------------
    def analyse(self):
        """Return (Report, {name: provenance}) from the data so far."""
        known = dict(self.given)
        prov = {k: "stated" for k in known}
        if self.model is not None and "T" in known and "P" in known:
            fluid, eos = self.model
            try:
                st = props.state(fluid, eos, known["T"], known["P"])
            except ThermoError as exc:
                self.say(f"  the EoS could not be evaluated: {exc}")
            else:
                for k in knowledge.QUANTITIES:
                    if k in st and k not in known:
                        known[k] = st[k]
                        prov[k] = f"parameterised ({self.model_name})"
        return knowledge.infer(known, prov), prov

    def pending(self, rep):
        """Targets that cannot be computed yet, with the reason."""
        out = []
        for t in self.targets:
            if t in self.skipped:
                continue
            try:
                if t in knowledge.QUANTITIES:
                    knowledge.require(rep, t)
                else:
                    knowledge.derivative_from_facts(t, rep)
            except ThermoError as exc:
                out.append((t, exc))
        return out

    def report(self, rep):
        if rep.violations or not rep.consistent:
            self.say("INCONSISTENT INPUT:")
            for v in rep.violations:
                self.say("  " + v)
            for rel, q, old, new, diff, ok in rep.checks:
                if not ok:
                    self.say(f"  {rel} gives {q} = {new:.6g}, but {q} = {old:.6g} was given "
                             f"(relative difference {diff:.3g})")
        self.say("Known and inferred quantities:")
        for name, f in rep.facts.items():
            self.say(f"  {name:7s} = {f.value:<14.6g} {SI_UNIT.get(name, ''):10s} [{f.provenance}]")
        self.say("Targets:")
        for t in self.targets:
            try:
                if t in knowledge.QUANTITIES:
                    v = knowledge.require(rep, t)
                    self.say(f"  {t} = {v:.6g} {SI_UNIT.get(t, '')}")
                else:
                    expr, v = knowledge.derivative_from_facts(t, rep)
                    self.say(f"  {t} = {expr} = {v:.6g}")
            except ThermoError as exc:
                self.say(f"  {t}: NOT COMPUTED. {exc}")

    # ---- the loop -------------------------------------------------------
    def run(self):
        self.say("thermocalculus interactive session. Type 'help' for commands.")
        self.say("First give the targets and the known values; the engine will ask for what is missing.")
        while True:
            rep, _ = self.analyse()
            todo = self.pending(rep) if self.targets and self.given else []
            if todo:
                target, exc = todo[0]
                self.say(f"\nCannot compute {target} yet: {exc}")
                prompt = "Enter a value (name = value unit), or 'skip' to give up on this target > "
            else:
                prompt = "> "
            line = self.read_line(prompt)
            if line is None:
                break
            line = line.strip()
            low = line.lower()
            if not line:
                continue
            if low in ("quit", "exit", "q"):
                break
            if low == "help":
                self.say(HELP)
            elif low == "show" or low == "go":
                self.report(rep)
            elif low == "skip" and todo:
                self.skipped.add(todo[0][0])
                self.say(f"  giving up on {todo[0][0]}")
            elif low.startswith("target"):
                for item in line[6:].split(";"):
                    self.add_target(item)
            elif low.startswith("model"):
                self.set_model(line[5:])
            elif "=" in line:
                self.add_value(line)
            else:
                self.say("  not understood; type 'help'")
        rep, _ = self.analyse()
        self.say("\nFinal result")
        self.report(rep)
        return rep



def run_session(script=None, fluid=None, eos=None, output_fn=print):
    """Start a session.  ``script`` is a path of pre-written answers (one per
    line, '#' starts a comment) used instead of the keyboard; each answer is
    echoed so that the transcript reads like the live session."""
    session = Session(output_fn=output_fn)
    if script:
        with open(script) as fh:
            answers = [l.rstrip("\n") for l in fh
                       if l.strip() and not l.lstrip().startswith("#")]
        it = iter(answers)

        def scripted(prompt):
            try:
                a = next(it)
            except StopIteration:
                raise EOFError
            output_fn(prompt + a)
            return a
        session.ask_user = scripted
    if fluid:
        session.set_model(f"{fluid} {eos or 'PR'}")
    return session.run()
