#!/usr/bin/env python3
"""
Prune a published, dated mammal/bat tree to the species you analysed and
relabel tips with your representative assembly names.

  prune_tree.py big_tree.nwk prdm9_qc.representatives.tsv -o bats
      [--map name_map.tsv]   species<TAB>tip_label overrides (synonyms!)

Tip matching: 'Genus species' matches a tip whose label starts with
Genus_species (case-insensitive), so Upham et al. 2019-style labels such as
Myotis_myotis_VESPERTILIONIDAE_CHIROPTERA work. Taxonomy has moved for
several bats in your table (e.g. Neoeptesicus, Cnephaeus, Dermanura,
Gardnerycteris, Lyroderma, Doryrhina) - unmatched species are listed in
<o>.missing.txt; add them to --map with the name the tree uses.

Outputs: <o>.species.nwk (tips = Genus_species), <o>.short.nwk (tips =
representative assembly short names, matching names written by
prdm9_downstream.sh for HyPhy), <o>.missing.txt
"""
import argparse, csv, re, sys


class Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self):
        self.name, self.length, self.children, self.parent = "", None, [], None


def parse(s):
    s = s.strip().rstrip(";")
    root = cur = Node()
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == "(":
            ch = Node(); ch.parent = cur; cur.children.append(ch); cur = ch; i += 1
        elif c == ",":
            ch = Node(); ch.parent = cur.parent; cur.parent.children.append(ch); cur = ch; i += 1
        elif c == ")":
            cur = cur.parent; i += 1
        elif c == ":":
            m = re.match(r":([^,();\[]+)", s[i:])
            cur.length = float(m.group(1)); i += len(m.group(0))
        elif c == "[":
            i = s.index("]", i) + 1
        elif c == "'":
            j = s.index("'", i + 1); cur.name = s[i + 1:j]; i = j + 1
        else:
            m = re.match(r"[^,();:\[]+", s[i:])
            cur.name = m.group(0).strip(); i += len(m.group(0))
    return root


def tips(node):
    if not node.children:
        return [node]
    return [t for c in node.children for t in tips(c)]


def prune(node, keep):
    if not node.children:
        return node if node in keep else None
    kids = [k for k in (prune(c, keep) for c in node.children) if k is not None]
    if not kids:
        return None
    if len(kids) == 1:                       # collapse unary node
        k = kids[0]
        if node.length is not None or k.length is not None:
            k.length = (k.length or 0) + (node.length or 0)
        return k
    node.children = kids
    for k in kids:
        k.parent = node
    return node


def write(node, label):
    if not node.children:
        s = label(node)
    else:
        s = "(" + ",".join(write(c, label) for c in node.children) + ")"
    if node.length is not None:
        s += f":{node.length:.6g}"
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tree")
    ap.add_argument("representatives")
    ap.add_argument("--map")
    ap.add_argument("-o", "--out", default="species_tree")
    a = ap.parse_args()

    text = open(a.tree).read()
    if "TREE" in text.upper() and "BEGIN" in text.upper():   # NEXUS: take first tree
        m = re.search(r"tree\s+[^=]+=\s*(\[[^\]]*\]\s*)?(\(.*?;)", text, re.I | re.S)
        text = m.group(2)
    root = parse(text)
    all_tips = tips(root)
    index = {}
    for t in all_tips:
        key = "_".join(t.name.replace(" ", "_").split("_")[:2]).lower()
        index.setdefault(key, t)
    by_label = {t.name: t for t in all_tips}

    override = {}
    if a.map:
        for l in open(a.map):
            if l.strip() and not l.startswith("#"):
                sp, lab = l.rstrip("\n").split("\t")[:2]
                override[sp] = lab
    reps = list(csv.DictReader(open(a.representatives), delimiter="\t"))
    keep, label_short, label_sp, missing = set(), {}, {}, []
    for r in reps:
        sp = r["species"]
        t = by_label.get(override[sp]) if sp in override else \
            index.get("_".join(sp.split()[:2]).lower())
        if t is None:
            missing.append(sp)
            continue
        if t in keep:
            print(f"note: {sp} maps to an already-used tip {t.name}", file=sys.stderr)
            continue
        keep.add(t)
        label_short[id(t)] = r["short"]
        label_sp[id(t)] = sp.replace(" ", "_")
    pr = prune(root, keep)
    pr.length = None                         # no root branch
    open(a.out + ".species.nwk", "w").write(write(pr, lambda t: label_sp[id(t)]) + ";\n")
    open(a.out + ".short.nwk", "w").write(write(pr, lambda t: label_short[id(t)]) + ";\n")
    open(a.out + ".missing.txt", "w").write("\n".join(missing) + ("\n" if missing else ""))
    print(f"kept {len(keep)} tips; {len(missing)} species unmatched -> {a.out}.missing.txt",
          file=sys.stderr)


if __name__ == "__main__":
    main()
