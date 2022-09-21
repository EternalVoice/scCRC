
import sys
import umap
import pandas as pd
import loompy as lp
from pyscenic.binarization import binarize
from pyscenic.cli.utils import load_signatures
from MulticoreTSNE import MulticoreTSNE as TSNE


f_pyscenic_output = sys.argv[1]
regulon_file = sys.argv[2]
project = sys.argv[3]
threads = int(sys.argv[4])
min_regulon_size = int(sys.argv[5])

def get_motif_logo(regulon):
    base_url = "http://motifcollections.aertslab.org/v9/logos/"
    for elem in regulon.context:
        if elem.endswith('.png'):
            return(base_url + elem)

# regulons to gmt file
print('''
##############################################
    1. Transform regulons to gmt file ...
##############################################
    ''')
regulons = load_signatures(regulon_file)
select_cols = [i.name for i in regulons if len(i.genes) >= min_regulon_size]
gmt_file = project + ".regulons.gmt"
txt_file = project + ".regulons.txt"
fo1 = open(gmt_file, 'w') 
fo2 = open(txt_file, 'w')
for i in regulons:
    if i.name in select_cols:
        motif = get_motif_logo(i)
        genes = "\t".join(i.genes)
        tf = "%s(%sg)" % (i.transcription_factor, len(i.genes))
        fo1.write("%s\t%s\t%s\n" % (tf, motif, genes))
        fo2.write("%s\t%s\t%s\n" % (tf, motif, genes.replace("\t",",")))

# collect SCENIC AUCell output
print('''
##############################################
    2. Collect SCENIC AUCell output ...
##############################################
    ''')
lf = lp.connect(f_pyscenic_output, mode='r+', validate=False)
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
auc_mtx = auc_mtx[select_cols]
auc_mtx.to_csv(project + ".AUCell.txt", sep='\t')
lf.close()

# Generate a binary regulon activity matrix
print('''
######################################################
    3. Generate a binary regulon activity matrix ...
######################################################
    ''')
binary_mtx, auc_thresholds = binarize(auc_mtx, num_workers=threads)
binary_mtx.to_csv(project + ".binary_mtx.txt", sep='\t')
auc_thresholds.to_csv(project + ".auc_thresholds.txt", sep='\t')

# UMAP
print('''
#######################
    4. UMAP ...
#######################
    ''')
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap(auc_mtx)
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv(project + ".umap.txt", sep='\t')
# tSNE
print('''
#######################
    5. tSNE ...
#######################
    ''')
tsne = TSNE( n_jobs=20 )
dr_tsne = tsne.fit_transform( auc_mtx )
pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv(project + ".tsne.txt", sep='\t')
