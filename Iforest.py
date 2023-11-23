import os
from typing import Literal
import pandas as pd
from sklearn.ensemble import IsolationForest
import re


def get_annotated(annotations, srna, strand=None):
    overlapping_annotations = annotations[
        annotations['start'].between(srna['start'], srna['end'], inclusive='neither') |
        annotations['end'].between(srna['start'], srna['end'], inclusive='neither') |
        ((annotations['start'] < srna['start']) & (srna['end'] < annotations['end']))
        ]
    if len(overlapping_annotations.index) == 0:
        return 'Intergenic'

    overlapped_rows = {}
    srna_length = srna['end'] - srna['start']
    for i, annotation in overlapping_annotations.iterrows():
        overlap = min(annotation['end'], srna['end']) - max(annotation['start'], srna['start'])
        overlapped_rows[i] = overlap * 100 / srna_length

    max_overlap = max(overlapped_rows.values())
    Is = [k for k, v in overlapped_rows.items() if v == max_overlap]

    if len(Is) == 1 or strand is None:
        I = Is[0]
    elif annotations.loc[Is[0], 'strand'] == strand:
        I = Is[0]
    elif annotations.loc[Is[1], 'strand'] == strand:
        I = Is[1]
    else:
        I = Is[0]

    if overlapped_rows[I] < 50 and len(Is) == 1:
        return 'Intergenic'

    match = re.fullmatch('^ID=([^;]*);.*$', annotations.loc[I, 'gene_id'])
    if match is None:
        return None
    elif strand is None:
        return match.group(1)
    elif annotations.loc[I, 'strand'] == strand:
        return 'Sense to ' + match.group(1)
    else:
        return 'Antisense to ' + match.group(1)


def annotate(regions, annotations, strand=None):
    d_out = pd.DataFrame(columns=['name', 'start', 'end', 'location'])
    for i, srna in regions.iterrows():
        an = get_annotated(annotations, srna, strand)
        name = 'sRNA-{0}'.format(i) if strand is None \
            else 'Fwd-sRNA-{0}'.format(i) if strand == '+' \
            else 'Rev-sRNA-{0}'.format(i)
        d_out.loc[i] = [name, srna['start'], srna['end'], an]
    return d_out


def group_consecutive(vals, step=1):
    run = []
    result = [run]
    expect = None
    for v in vals:
        if (v == expect) or (expect is None):
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step
    return result


def detectSRNA(data, cf):
    readcount = data['readcount'].values.reshape(-1, 1)
    clf = IsolationForest(contamination=cf, n_jobs=1, n_estimators=100, max_samples=5000, random_state=0)
    clf.fit(readcount)

    pred = clf.predict(readcount)
    scores = clf.decision_function(readcount)

    # add the anomaly score and decision(-1/1) to the data frame
    data['anomaly decision'] = pred
    data['score'] = scores

    # Get the regions
    a = data[data['anomaly decision'] == -1]
    Regions = []
    for group in group_consecutive(a['position']):
        if len(group) > 1:
            Regions.append([group[0], group[-1]])

    Merged = [Regions[0]]
    for region in Regions[1:]:
        if region[0] - Merged[-1][-1] <= 5:
            Merged[-1][-1] = region[-1]
        else:
            Merged.append(region)

    reg = pd.DataFrame(Merged, columns=['start', 'end'])
    regions = reg[reg['end'] - reg['start'] >= 20].copy()
    return regions


org_names = Literal[
    'Bacillus subtilis 168',
    'Clostridium difficile S-0253',
    'Escherichia coli K-12',
    'Listeria monocytogenes EGD-e',
    'Mycobacterium tuberculosis H37Rv',
    'Pseudomonas aeruginosa PAO1',
    'Salmonella enterica LT2',
    'Staphylococcus aureus NCTC 8325',
    'Streptococcus pneumoniae Hu17',
    'Vibrio cholerae RFB16']
strand_types = Literal[
    'strand_spec',
    'non_strand_spec'
]


def run_posif(
        organism_name: org_names,
        strand_type: strand_types,
        output_dir: str,
        fwd_perbase_file: str = None,
        rev_perbase_file: str = None,
        perbase_file: str = None,
        contamination_factor: float = 0.05,
):
    anno_file = os.path.join('annotations', '{0}.gff'.format(organism_name))
    annotation_data = pd.read_csv(anno_file, sep='\t', comment='#', usecols=[2, 3, 4, 6, 8],
                                  names=['type', 'start', 'end', 'strand', 'gene_id'])
    gene_annotations = annotation_data[annotation_data['type'] == 'gene'].copy()
    del annotation_data

    if strand_type == 'strand_spec':
        if fwd_perbase_file is None or rev_perbase_file is None:
            raise ValueError('For Strand-Specific, fwd_perbase_file as well as '
                             'rev_perbase_file should be provided')
        fwd_data = pd.read_csv(fwd_perbase_file, sep='\t', names=['position', 'readcount'], usecols=[1, 2])
        rev_data = pd.read_csv(rev_perbase_file, sep='\t', names=['position', 'readcount'], usecols=[1, 2])
        fwd_regions = detectSRNA(fwd_data, contamination_factor)
        rev_regions = detectSRNA(rev_data, contamination_factor)
        fwd_output = annotate(fwd_regions, gene_annotations, '+')
        rev_output = annotate(rev_regions, gene_annotations, '-')

        fwd_output.to_csv(os.path.join(output_dir, 'sRNA_fwd_output.csv'), index=False)
        rev_output.to_csv(os.path.join(output_dir, 'sRNA_rev_output.csv'), index=False)

    else:
        if perbase_file is None:
            raise ValueError('For non Strand-Specific, perbase_file should be provided')
        perbase_data = pd.read_csv(perbase_file, sep='\t', names=['position', 'readcount'], usecols=[1, 2])

        regions = detectSRNA(perbase_data, contamination_factor)
        output = annotate(regions, gene_annotations)
        output.to_csv(os.path.join(output_dir, 'sRNA_output.csv'), index=False)


# ======================= Example Run =========================
run_posif(
    organism_name='Mycobacterium tuberculosis H37Rv',
    strand_type='strand_spec',
    fwd_perbase_file='Exampleinput/Control1F.bed',
    rev_perbase_file='Exampleinput/Control1R.bed',
    output_dir='Exampleoutput',
    contamination_factor=0.05
)
# =============================================================
