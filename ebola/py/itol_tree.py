from pastml import numeric2datetime
from pastml.tree import read_forest, DATE


HEADER = """DATASET_COLORSTRIP

SEPARATOR TAB
BORDER_WIDTH	0.5
MARGIN	5
STRIP_WIDTH	25
COLOR	#ffffff
DATASET_LABEL	in_dataset

LEGEND_COLORS	#d53e4f
LEGEND_LABELS	selected
LEGEND_SHAPES	1
LEGEND_TITLE	in_dataset

DATA
#NODE_ID COLOR LABEL_OR_STYLE
"""
if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--tree', required=True, type=str)
    parser.add_argument('--itol', required=True, type=str)
    params = parser.parse_args()

    forest = read_forest(params.tree)
    with open(params.itol, 'w+') as f:
        f.write(HEADER)
        for tree in forest:
            for t in tree:
                f.write('{}\t#d53e4f\tselected\n'.format(t.name, ))
