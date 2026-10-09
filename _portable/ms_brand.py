"""Presentation names; internal engine identifiers and licences remain intact."""
import re

# (pattern, replacement), in this order: engine phrases first, so "the UniDec
# engine" does not become "the Bayesian deconvolution engine"; "Bayesian
# engine" is the name 3.4 to 3.55 gave (saved results and projects keep it).
# Applied to what the windows show: method, notes, errors, progress, quality.
_NAMES = [(re.compile(p), r) for p, r in (
    (r"UniDec's IsoDec", "IsoDec"),
    (r"UniDec \(Bayesian\)", "Bayesian deconvolution"),
    (r"MetaUniDec", "batch data"),
    (r"(^|[.:;!?]\s+)UniDec's (engine|scores)", r"\1The \2"),
    (r"UniDec's (engine|scores)", r"the \1"),
    (r"([Tt]he )?UniDec(?: [\d.]+)? engine", r"\1deconvolution engine"),
    (r"the UniDec default", "the default"),
    (r"UniDec cannot deconvolute this spectrum", "This spectrum cannot be deconvoluted"),
    (r"UniDec", "Bayesian deconvolution"),
    (r"([Tt]he )?Bayesian engine engine", r"\1deconvolution engine"),
    (r"Bayesian engine", "Bayesian deconvolution"),
)]


def display(text):
    """The text with the engine named as in the program's windows."""
    if not isinstance(text, str):
        return text
    for pat, rep in _NAMES:
        text = pat.sub(rep, text)
    return text


CREDITS=('Third-party software: MSpektra uses UniDec. Please cite Marty et al., '
         'Analytical Chemistry 2015, DOI: 10.1021/acs.analchem.5b00140 when publishing results obtained with that engine. '
         'Copyright and licence notices are retained in the LICENSES folder.')
