from openpyxl import Workbook
from openpyxl.styles import Alignment


def parse(file_name):
    proteins = []
    di1 = {}
    with open(file_name) as fil:
        for line in fil.readlines():

            # name and size of protein
            if line.startswith('ID   '):
                if di1:
                    proteins.append(di1)
                    di1 = {}
                # print(line)
                _, name, _, length, _ = line.split()
                di1['name'] = name
                di1['size'] = length
                # print(name, length)

            # entry number
            elif line.startswith('AC   '):
                _, entry, *_ = line.split()
                # print(entry[:-1])
                di1['entry'] = entry[:-1]

            # protein name
            elif line.startswith('DE   RecName: Full='):
                # print(line[19:-2])
                di1['fullname'] = line[19:-2]

            # organism data
            elif line.startswith('OS   '):
                _, *organism = line.split()
                # print(' '.join(organism[:-1]))
                di1['organism'] = ' '.join(organism[:-1])

            # transmembrane data
            elif line.startswith('FT   TRANSMEM'):
                _, _, rng = line.split()
                apo, eos = rng.split('..')
                # print(name, apo, eos)
                apo, eos = int(apo.replace('>', '').replace('<', '')), int(
                    eos.replace('>', '').replace('<', ''))
                di1['hreg'] = di1.get('hreg', [])
                di1['hreg'].append({'apo': apo, 'eos': eos})

            elif line.startswith('DR   PDB;'):
                # print(line)
                _, tag, method, version, rest = line.split(';')
                tag = tag.strip()
                method = method.strip()
                rest = rest.strip()[:-1]

                # print(method, version, rest)
                if method == 'Model':
                    continue
                if not rest:
                    continue
                rregion = rest.split(',')

                # print(rregion)
                for region in rregion:
                    rcategory, apoeos = region.split('=')
                    rapo, reos = apoeos.split('-')
                    rcategory = rcategory.strip()
                    rapo = int(rapo.strip())
                    reos = int(reos.strip())
                    di1['cref'] = di1.get('cref', [])
                    di1['cref'].append({'tag': tag, 'method': method, 'ver': version,
                                        'cat': rcategory, 'apo': int(rapo), 'eos': int(reos)})

    proteins.append(di1)
    return proteins


def contains(cstart, cend, sstart, send):
    if cstart <= sstart and cend >= send:
        return True
    return False


def dic_contains(dic1):
    if 'cref' not in dic1:
        dic1['error'] = 'No cross reference data for protein'
        dic1['cref'] = []
    if 'hreg' not in dic1:
        dic1['error'] = 'No reference to helical regions'
        dic1['hreg'] = []
    for trn in dic1['hreg']:
        trn['result'] = False
        for region in dic1['cref']:
            if contains(region['apo'], region['eos'], trn['apo'], trn['eos']):
                trn['result'] = True
                # print(trn, region)
                break


def print_results(filename):
    proteins = parse(filename)
    print('┌----------------------------------------------------------------------------------------┐')
    print('| entry\t\t', ' prot_name\t\t',
          'cover_status\t\t', 'cross_referenced_regions|')
    print('└----------------------------------------------------------------------------------------┘')
    for protein in proteins:
        dic_contains(protein)
        lhreg = len(protein['hreg'])
        lcref = sum([p['result'] for p in protein['hreg']])
        protein['number_of_transmembrane_regions'] = lhreg
        protein['number_of_cross_referenced_regions'] = lcref
        if lcref == 0:
            protein['cover_status'] = 'none'
        elif lhreg == lcref:
            protein['cover_status'] = 'all'
        else:
            protein['cover_status'] = 'some'
        print(
            f"| {protein['entry']:15} {protein['name']:26} {protein['cover_status']:32} {protein['number_of_cross_referenced_regions']}\t\t |")
    print('└----------------------------------------------------------------------------------------┘')


def write_xl(protein_file, xlfile):
    proteins = parse(protein_file)
    wb = Workbook()
    sheet = wb.active
    sheet["A1"] = "Entry"
    sheet["B1"] = "Protein name"
    sheet["C1"] = "Cover status"
    sheet["D1"] = "Number of cross referenced regions"
    sheet["A1"].alignment = Alignment(
        wrap_text=True, vertical='center', horizontal='center')
    sheet["B1"].alignment = Alignment(
        wrap_text=True, vertical='center', horizontal='center')
    sheet["C1"].alignment = Alignment(
        wrap_text=True, vertical='center', horizontal='center')
    sheet["D1"].alignment = Alignment(
        wrap_text=True, vertical='center', horizontal='center')
    sheet.row_dimensions[1].height = 70
    sheet.column_dimensions['A'].width = 12
    sheet.column_dimensions['B'].width = 18
    sheet.column_dimensions['C'].width = 8
    sheet.column_dimensions['D'].width = 9
    for i, protein in enumerate(proteins):
        dic_contains(protein)
        lhreg = len(protein['hreg'])
        lcref = sum([p['result'] for p in protein['hreg']])
        protein['number_of_transmembrane_regions'] = lhreg
        protein['number_of_cross_referenced_regions'] = lcref
        if lcref == 0:
            protein['cover_status'] = 'none'
        elif lhreg == lcref:
            protein['cover_status'] = 'all'
        else:
            protein['cover_status'] = 'some'
        sheet[f"A{i+2}"] = protein['entry']
        sheet[f"B{i+2}"] = protein['name']
        sheet[f"C{i+2}"] = protein['cover_status']
        sheet[f"D{i+2}"] = protein['number_of_cross_referenced_regions']
    wb.save(filename=xlfile)


if __name__ == '__main__':
    print_results('g_proteins.txt')
    write_xl('g_proteins.txt', 'g_protein_results.xlsx')
