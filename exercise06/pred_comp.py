H = 'H'  # aplha-helix
E = 'E'  # beta-fold
T = 'C'  # coil
K = '-'  # uknown

STC = {'H': H, 'E': E, 'T': T, 'K': K}

def make_sec_sequence(data):
    """Create secondary sequence from uniprot secondary structure data"""
    empty = '*'
    final = ''
    pointer = 0
    for elm in data:
        start, end, type = elm
        size = end - start + 1
        stripe = size * type
        if pointer > start - 1:
            raise ValueError('Error')
        if pointer == (start - 1):
            final += stripe
        else:
            gap_size = start - 1 - pointer
            final += gap_size * K
            final += stripe
        pointer = end
    return final

def comparator(real, pred, typ='H'):
    """Compares real secondary structure data to predicted"""

    assert len(real) == len(pred)
    tpv = {'TruePositive': 0, 'TrueNegative': 0 , 'FalsePositive': 0, 'FalseNegative': 0, 'type': typ, 'q': 0}
    for i, realv in enumerate(real):
        if realv == typ:
            if pred[i] == typ:
                tpv['TruePositive'] += 1
            else:
                tpv['FalseNegative'] += 1
        else:
            if pred[i] == typ:
                tpv['FalsePositive'] += 1
            else:
                tpv['TrueNegative'] += 1
    tpv['q'] = cq1(tpv)
    return tpv


def cq1(vls):
    """Calculate q value"""
    nume = vls['TruePositive'] + vls['TrueNegative']
    denu = vls['TruePositive'] + vls['TrueNegative'] + vls['FalsePositive'] + vls['FalseNegative']
    return  round(nume / denu  * 100, 2)


def parse_horiz(fname):
    """Parse file from psi-pred horiz file type"""
    res = ''
    with open(fname) as fil:
        for lin in fil.readlines():
            if lin.startswith('Pred'):
                res += lin.strip()[6:]
    return res


def parse_jpred(fname):
    """Parse file from jpred prediction"""
    res = ''
    with open(fname) as fil:
        for lin in fil.readlines():
            if lin.startswith('PRED'):
                res += lin.strip()[6:]
    return res


def parse_porter(fname):
    """Parse file from porter prediction"""
    res = ''
    with open(fname) as fil:
        for lin in fil.readlines():
            if lin.startswith('PRED'):
                res += lin.strip()[6:]
    return res


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


def calc_q3(rel, pred):
    """Calculate q3 value"""
    res = []
    q3 = 0
    for i in 'HEC':
        res.append(comparator(rel, pred, i))
    for val in res:
        q3 += val['q']
    # print(res)
    return round(q3 / len(res), 2), res


def calc_q3_jprd_porter_psi(datdir, verbose=False):
    """Calculate q3 values from every program and print"""
    linesize = 50
    real = make_sec_sequence(parse_real(f'{datdir}/real.txt'))
    # print(real)
    jprd = parse_jpred(f'{datdir}/jpred.txt')
    porr = parse_porter(f'{datdir}/porter.txt')
    pspr = parse_horiz(f'{datdir}/psipred.horiz')
    q3_jpred, dt_jpred = calc_q3(real, jprd)
    q3_portr, dt_portr = calc_q3(real, porr)
    q3_psipr, dt_psipr = calc_q3(real, pspr)
    print(f'Q3 factor results per algorithm for {datdir} protein')
    print(f"Q3(jpred)    : {q3_jpred}")
    print(f"Q3(porter)   : {q3_portr}")
    print(f"Q3(psi-pred) : {q3_psipr}")
    if not verbose:
        return

    print('=' * linesize)

    print(f'\nAnalysis for protein: {datdir} and algorithm: jpred\n')
    for elm in dt_jpred:
        print(f"Secondary structure type: {elm['type']}  q: {elm['q']}")
        print(f"True positive : {elm['TruePositive']:4.0f}     False positive: {elm['FalsePositive']:4.0f}")
        print(f"False negative: {elm['FalseNegative']:4.0f}      True negative: {elm['TrueNegative']:4.0f}\n")

    print('-' * linesize)
    print(f'\nAnalysis for protein: {datdir} and algorithm: porter\n')
    for elm in dt_portr:
        print(f"Secondary structure type: {elm['type']}  q: {elm['q']}")
        print(f"True positive : {elm['TruePositive']:4.0f}     False positive: {elm['FalsePositive']:4.0f}")
        print(f"False negative: {elm['FalseNegative']:4.0f}      True negative: {elm['TrueNegative']:4.0f}\n")

    print('-' * linesize)
    print(f'\nAnalysis for protein: {datdir} and algorithm: psi-pred\n')
    for elm in dt_psipr:
        print(f"Secondary structure type: {elm['type']}  q: {elm['q']}")
        print(f"True positive : {elm['TruePositive']:4.0f}     False positive: {elm['FalsePositive']:4.0f}")
        print(f"False negative: {elm['FalseNegative']:4.0f}      True negative: {elm['TrueNegative']:4.0f}\n")

    print('=' * linesize)
    # print(dt_jpred)
    # print(dt_portr)
    # print(dt_psipr)


if __name__ == '__main__':
    calc_q3_jprd_porter_psi('ompa', True)
    print("")
    calc_q3_jprd_porter_psi('malic', True)

