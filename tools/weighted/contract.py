"""Contract low-support branches of gene trees into polytomies.

usage: python3 -I tools/weighted/contract.py -t THRESHOLD IN OUT [IN OUT ...]
       python3 -I tools/weighted/contract.py -t THRESHOLD --tree DATA_IN DATA_OUT [--conditions LIST] [--jobs N]

An internal branch whose support (the label after its ')') is below THRESHOLD is removed and its
children attach to its parent; branches without a numeric support are kept. Leaf labels, kept
supports and branch lengths are written unchanged, one tree per line. With --tree, every file named
estimatedgenetre under DATA_IN/estimated-gene-trees is contracted into the same path under
DATA_OUT (the run_all.py layout), and DATA_OUT/true_species_trees links to DATA_IN's.
"""

import argparse
import os
import re
import sys
from multiprocessing import Pool

sys.setrecursionlimit(100000)

TOKEN = re.compile(r'\(|\)|,|;|[^(),;]+')


def parse(text):
    """Nested lists: an internal node is [children, label]; a leaf is its label text (with length)."""
    stack = [[[], '']]
    prev = ''
    for tok in TOKEN.findall(text.strip()):
        if tok == '(':
            stack.append([[], ''])
        elif tok == ')':
            node = stack.pop()
            stack[-1][0].append(node)
        elif tok in ',;':
            pass
        elif prev == ')':
            stack[-1][0][-1][1] = tok           # the label (support[:length]) of the node just closed
        else:
            stack[-1][0].append(tok)            # a leaf: label[:length]
        prev = tok
    root = stack[0][0]
    if len(stack) != 1 or len(root) != 1:
        raise ValueError('unbalanced parentheses')
    return root[0]


def support(label):
    s = label.split(':', 1)[0]
    try:
        return float(s)
    except ValueError:
        return None


def contract(node, threshold):
    """Returns the list of nodes that replace `node` in its parent (itself, or its children)."""
    if isinstance(node, str):
        return [node]
    kids = []
    for c in node[0]:
        kids.extend(contract(c, threshold))
    node[0] = kids
    s = support(node[1])
    if s is not None and s < threshold:
        return kids
    return [node]


def write(node):
    if isinstance(node, str):
        return node
    return '(' + ','.join(write(c) for c in node[0]) + ')' + node[1]


def contract_file(args):
    src, dst, threshold = args
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    n = removed = 0
    with open(src) as f, open(dst + '.part', 'w') as out:
        for line in f:
            if '(' not in line:
                continue
            root = parse(line)
            before = line.count(')')
            root[0] = [x for c in root[0] for x in contract(c, threshold)]   # the root itself stays
            text = write(root)
            removed += before - text.count(')')
            out.write(text + ';\n')
            n += 1
    os.replace(dst + '.part', dst)
    return src, n, removed


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('-t', '--threshold', type=float, required=True)
    ap.add_argument('--tree', action='store_true', help='IN and OUT are data folders in the run_all.py layout')
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--conditions', help='with --tree: only these conditions (comma-separated)')
    ap.add_argument('paths', nargs='+')
    a = ap.parse_args()
    if len(a.paths) % 2:
        sys.exit('paths come in IN OUT pairs')
    jobs = []
    for src, dst in zip(a.paths[::2], a.paths[1::2]):
        if a.tree:
            gi = os.path.join(src, 'estimated-gene-trees')
            conds = a.conditions.split(',') if a.conditions else sorted(os.listdir(gi))
            for cond in conds:
                for rep in sorted(os.listdir(os.path.join(gi, cond))):
                    f = os.path.join(gi, cond, rep, 'estimatedgenetre')
                    if os.path.isfile(f):
                        jobs.append((f, os.path.join(dst, 'estimated-gene-trees', cond, rep, 'estimatedgenetre'), a.threshold))
            os.makedirs(dst, exist_ok=True)
            link = os.path.join(dst, 'true_species_trees')
            if not os.path.exists(link):
                os.symlink(os.path.realpath(os.path.join(src, 'true_species_trees')), link)
        else:
            jobs.append((src, dst, a.threshold))
    total = removed = 0
    with Pool(a.jobs) as pool:
        for src, n, r in pool.imap_unordered(contract_file, jobs):
            total += n
            removed += r
    print(f'threshold {a.threshold}: {len(jobs)} files, {total} trees, {removed} branches contracted')


if __name__ == '__main__':
    main()
