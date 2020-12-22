I = 'I'  # cytoplasmic (inside)
M = 'M'  # transmembrane
O = 'O'  # extracellular (outside)
K = '-'  # uknown

STC = {'I': I, 'M': M, 'O': O, 'K': K}


def make_trans_sequence(data):
    """Create secondary sequence from uniprot secondary structure data"""
    final = ''
    pointer = 0
    for elm in data:
        start, end, type = elm
        size = end - start + 1
        stripe = size * type
        if pointer > start - 1:
            raise ValueError(f'Pointer({pointer}) > {start - 1}')
        if pointer == (start - 1):
            final += stripe
        else:
            gap_size = start - 1 - pointer
            final += gap_size * K
            final += stripe
        pointer = end
    return final


# def parse_horiz(fname):
#     """Parse file from psi-pred horiz file type"""
#     res = ''
#     with open(fname) as fil:
#         for lin in fil.readlines():
#             if lin.startswith('Pred'):
#                 res += lin.strip()[6:]
#     return res


# def parse_jpred(fname):
#     """Parse file from jpred prediction"""
#     res = ''
#     with open(fname) as fil:
#         for lin in fil.readlines():
#             if lin.startswith('PRED'):
#                 res += lin.strip()[6:]
#     return res


# def parse_porter(fname):
#     """Parse file from porter prediction"""
#     res = ''
#     with open(fname) as fil:
#         for lin in fil.readlines():
#             if lin.startswith('PRED'):
#                 res += lin.strip()[6:]
#     return res


def parse_real(fname):
    """Parse file with uniprot/validated real data"""
    res = []
    with open(fname) as fil:
        for lin in fil.readlines():
            slin = lin.strip()
            if len(slin) < 5:
                continue
            apo, eos, typ = slin.split()
            res.append((int(apo), int(eos), STC[typ]))
    return res


if __name__ == '__main__':
    datdir = 'rhodopsin'
    real = make_trans_sequence(parse_real(f'{datdir}/real.txt'))
    thmm = make_trans_sequence(parse_real(f'{datdir}/thmm.txt'))
    memsat = make_trans_sequence(parse_real(f'{datdir}/memsat.txt'))
    phobius = make_trans_sequence(parse_real(f'{datdir}/phobius.txt'))
    hmmtop = make_trans_sequence(parse_real(f'{datdir}/hmmtop.txt'))
    topcons = make_trans_sequence(parse_real(f'{datdir}/topcons.txt'))
    # print("")
    # print(real)
    # print("")
    print(thmm)
    print("")
    print(memsat)
    print("")
    print(phobius)
    print("")
    print(hmmtop)
    print("")
    print(topcons)
