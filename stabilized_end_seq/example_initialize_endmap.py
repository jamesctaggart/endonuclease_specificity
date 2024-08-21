from endmap.datastructures import EndMap

mochi_file = 'CDS_168.mochiview_v13.txt'

from endmap.datastructures import EndMap

# Initialize an EndMap object from .wig files.
def load_EM(wigdir, wigf, wigr, norm=True):
    print(wigf)
    out = EndMap()
    out.read_wigs(wigdir, wigf, wigr)
    if norm:
        print('RPM...')
        out.calculate_RPM()
        print('CDS...')
        out.normalize_to_cds(mochi_file)  # This is an annotation file specifying positions fo genes

    return out

em = load_EM(datadir, 'wig_5f_filename', 'wig_5r_filename') # 5f = 5' mapped, forward strand