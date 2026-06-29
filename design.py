"""
GFP thermostable variant designer for Synbio Challenges 2026 Protein Design track.
Reproducible pipeline for BioSyn-GP-Team.

Usage:
    python design.py

Output:
    submission.csv      - 6 designed sequences
    design_log.csv      - per-sequence design log
"""

import os
import itertools
import pandas as pd
import numpy as np
import torch

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

REFERENCE_FASTA = os.path.join(DATA_DIR, 'AAseqs of 5 GFP proteins_20260511.txt')
EXCLUSION_CSV = os.path.join(DATA_DIR, 'Exclusion_List.csv')

TEAM_NAME = 'BioSyn-GP-Team'
TARGET_LEN_MIN = 220
TARGET_LEN_MAX = 250


def load_reference_sequences(path):
    seqs = {}
    cur = None
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                cur = line[1:].split()[0]
                seqs[cur] = ''
            elif line and not line.startswith('#'):
                seqs[cur] += line
    return seqs


def load_exclusion_set(path):
    df = pd.read_csv(path)
    return set(df['Sequence'])


def mutate(seq, mutations):
    """Apply list of (pos, new_aa) 1-indexed mutations to seq."""
    s = list(seq)
    for pos, aa in mutations:
        s[pos - 1] = aa
    return ''.join(s)


def build_candidates(backbone, mutation_pool, excl_set):
    """Generate all combinatorial variants of size 3-7 from mutation_pool."""
    candidates = []
    for k in range(3, min(8, len(mutation_pool) + 1)):
        for combo in itertools.combinations(mutation_pool, k):
            var = mutate(backbone, combo)
            if TARGET_LEN_MIN <= len(var) <= TARGET_LEN_MAX:
                if var not in excl_set and var != backbone:
                    candidates.append((list(combo), var))
    return candidates


def score_mutations(muts, beneficial=None):
    """Score variant by mutation count + beneficial-site bonus."""
    if beneficial is None:
        beneficial = {(65, 'T'), (69, 'L'), (148, 'D'), (163, 'A'),
                      (203, 'F'), (222, 'L'), (231, 'L')}
    n = len(muts)
    bonus = sum(1 for m in muts if m in beneficial)
    penalty = max(0, (n - 7)) * 2
    return n * 2 + bonus - penalty


def compute_esm_embedding_norm(sequences, model, alphabet, batch_converter):
    """Compute mean embedding L2 norm per sequence via ESM-2."""
    data = [(f'seq_{i}', s) for i, s in enumerate(sequences)]
    _, _, batch_tokens = batch_converter(data)
    device = next(model.parameters()).device
    batch_tokens = batch_tokens.to(device)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[6])
    embeddings = results['representations'][6][:, 1:-1].mean(dim=1)
    return embeddings.norm(dim=1).cpu().numpy()


def main():
    print('=' * 60)
    print('BioSyn-GP-Team · GFP Thermostable Variant Designer')
    print('=' * 60)

    # --- Step 1: Load data ---
    print('\n[Step 1] Loading reference sequences...')
    seqs = load_reference_sequences(REFERENCE_FASTA)
    for name, seq in seqs.items():
        print(f'  {name}: len={len(seq)}')
    sfGFP = seqs['sfGFP']
    avGFP = seqs['avGFP']
    amacGFP = seqs['amacGFP']

    print('\n[Step 2] Loading exclusion list...')
    excl_set = load_exclusion_set(EXCLUSION_CSV)
    print(f'  Size: {len(excl_set)}')

    # --- Step 2: Define mutation pools ---
    print('\n[Step 3] Defining mutation pools...')
    avGFP_full = [(65, 'T'), (69, 'L'), (148, 'D'), (163, 'A'),
                  (198, 'S'), (203, 'F'), (222, 'L'), (231, 'L')]
    sfGFP_full = [(65, 'T'), (69, 'L'), (148, 'D'), (163, 'A'),
                  (198, 'S'), (203, 'F'), (222, 'L'), (231, 'L'),
                  (2, 'C'), (39, 'Y')]
    amacGFP_full = [(65, 'T'), (148, 'D'), (163, 'A'), (198, 'S'),
                    (203, 'F'), (222, 'L'), (231, 'L')]

    avGFP_pool = [(p, a) for (p, a) in avGFP_full if avGFP[p-1] != a]
    sfGFP_pool = [(p, a) for (p, a) in sfGFP_full if sfGFP[p-1] != a]
    amacGFP_pool = [(p, a) for (p, a) in amacGFP_full if amacGFP[p-1] != a]
    print(f'  avGFP pool ({len(avGFP_pool)}): {avGFP_pool}')
    print(f'  sfGFP pool ({len(sfGFP_pool)}): {sfGFP_pool}')
    print(f'  amacGFP pool ({len(amacGFP_pool)}): {amacGFP_pool}')

    # --- Step 3: Generate and filter candidates ---
    print('\n[Step 4] Generating combinatorial candidates...')
    variants = []
    for backbone, pool, name in [
        (avGFP, avGFP_pool, 'avGFP'),
        (sfGFP, sfGFP_pool, 'sfGFP'),
        (amacGFP, amacGFP_pool, 'amacGFP'),
    ]:
        cands = build_candidates(backbone, pool, excl_set)
        for muts, var in cands:
            variants.append({
                'name': f'{name}+' + '+'.join(f'{p}{a}' for p, a in muts),
                'seq': var,
                'muts': muts,
                'backbone': name,
            })
        print(f'  {name}: {len(cands)} candidates after filtering')

    # Deduplicate
    seen = set()
    unique = []
    for v in variants:
        if v['seq'] not in seen:
            seen.add(v['seq'])
            unique.append(v)
    variants = unique
    print(f'  Total unique candidates: {len(variants)}')

    # Score by mutation count/quality
    for v in variants:
        v['mut_score'] = score_mutations(v['muts'])
    variants.sort(key=lambda v: v['mut_score'], reverse=True)

    # --- Step 4: Diversity quota selection ---
    print('\n[Step 5] Diversity-based selection...')
    selected = []
    per_backbone_limit = 4
    counts = {'avGFP': 0, 'sfGFP': 0, 'amacGFP': 0}
    for v in variants:
        b = v['backbone']
        if counts[b] < per_backbone_limit:
            selected.append(v)
            counts[b] += 1
        if len(selected) >= 6:
            break
    for v in selected:
        print(f'  {v["name"]}: mut_score={v["mut_score"]}, len={len(v["seq"])}')

    # --- Step 5: ESM-2 validation ---
    print('\n[Step 6] ESM-2 (8M) embedding validation...')
    try:
        import esm
        model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
        model.eval()
        batch_converter = alphabet.get_batch_converter()
        seqs_to_score = [v['seq'] for v in selected]
        norms = compute_esm_embedding_norm(seqs_to_score, model, alphabet, batch_converter)
        for v, n in zip(selected, norms):
            v['esm_emb_norm'] = float(n)
            print(f'  {v["name"]}: esm_emb_norm={n:.3f}')
        selected.sort(key=lambda v: v['esm_emb_norm'], reverse=True)
        esm_valid = True
    except Exception as e:
        print(f'  ESM-2 validation skipped: {e}')
        for v in selected:
            v['esm_emb_norm'] = None
        esm_valid = False

    final = selected[:6]

    # --- Step 6: Final validation ---
    print('\n[Step 7] Final validation...')
    for v in final:
        seq = v['seq']
        assert seq[0] == 'M', f'{v["name"]}: must start with M'
        assert all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in seq), f'{v["name"]}: invalid AA'
        assert TARGET_LEN_MIN <= len(seq) <= TARGET_LEN_MAX, f'{v["name"]}: length'
        assert seq not in excl_set, f'{v["name"]}: in exclusion list!'
        print(f'  ✓ {v["name"]}: len={len(seq)}')

    # --- Step 7: Write outputs ---
    print('\n[Step 8] Writing outputs...')
    rows = [{
        'Team_Name': TEAM_NAME,
        'Seq_ID': i + 1,
        'Sequence': v['seq'],
    } for i, v in enumerate(final)]
    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, 'submission.csv'), index=False)

    log_rows = [{
        'Seq_ID': i + 1,
        'Name': v['name'],
        'Backbone': v['backbone'],
        'Mutations': str(v['muts']),
        'Mutation_Score': v['mut_score'],
        'ESM2_Embedding_Norm': v.get('esm_emb_norm'),
        'Length': len(v['seq']),
    } for i, v in enumerate(final)]
    pd.DataFrame(log_rows).to_csv(os.path.join(OUT_DIR, 'design_log.csv'), index=False)

    print(f'  ✓ submission.csv saved')
    print(f'  ✓ design_log.csv saved')
    print('\n' + '=' * 60)
    print('DONE. 6 sequences ready for submission.')
    print('=' * 60)


if __name__ == '__main__':
    main()
