"""
Sign rules on the ecological Grassmannian: decide impossibility of a
Lotka–Volterra sign pattern by propagating signs of the minors of B.

Typical use::

    from signrules import Ecology, complete

    eco = Ecology.parse("+---", "+-++/-+++/+++-/++-+")   # signs of a; rows of B (b_ij < 0: j benefits i)
    print(eco)                                           # the sign matrix [diag(a) | B]
    result = complete(eco)
    print(result.verdict, result.summary())              # impossible | open | undecided, and one line why

    from signrules import ELEVEN                         # the named patterns: ELEVEN, PATTERNS, presets()
    a, B = ELEVEN["S-a obligate mutualism"]
    print(complete(Ecology(a, B)).summary())             # "impossible at the root (R3 balanced ...)"

See docs/SIGNRULES.md for the mathematics and rules.py for the rule catalogue.
"""
from .plucker_index import PluckerIndex, UNKNOWN, index
from .chirotope import build_relations, build_relation_ids, propagate, violates
from .rules import Ecology, Conflict
from .search import complete, SearchResult, Derivation
from .certificate import (certificate_diagram, certificate_dot, write_certificate,
                          completion_diagram, completion_dot, write_completion,
                          search_tree_text, search_tree_diagram, search_tree_dot, write_search_tree)
from .gvrender import Diagram
from . import labels, style
from .treeplot import plot_search_tree, save_search_tree_png, search_tree_plotly, save_search_tree_html
from .replay import Replayer, ReplayReport, replay_certificate, replay_kill
from . import patterns
from .patterns import (ELEVEN, PATTERNS, parity_pattern, random_pattern, canonical_classes, trivial,
                       parse_signs, parse_pattern, presets, resolve_preset)

__all__ = ["PluckerIndex", "UNKNOWN", "index", "build_relations", "build_relation_ids", "propagate", "violates",
           "Ecology", "Conflict", "complete", "SearchResult", "Derivation",
           "Replayer", "ReplayReport", "replay_certificate", "replay_kill",
           "certificate_diagram", "certificate_dot", "write_certificate",
           "completion_diagram", "completion_dot", "write_completion",
           "search_tree_text", "search_tree_diagram", "search_tree_dot", "write_search_tree",
           "Diagram", "labels", "style", "patterns",
           "ELEVEN", "PATTERNS", "parity_pattern", "random_pattern", "canonical_classes", "trivial",
           "parse_signs", "parse_pattern", "presets", "resolve_preset",
           "plot_search_tree", "save_search_tree_png", "search_tree_plotly", "save_search_tree_html"]
