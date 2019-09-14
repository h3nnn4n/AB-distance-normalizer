import requests
import gzip
import sys
import csv
import os


def fetch_native_pdb(target):
    if not os.path.isfile(target + '.pdb.gz'):
        print('[INFO] Fetching %s from PDB' % target)
        url = "http://files.rcsb.org/download/%s.pdb.gz" % target
        r = requests.get(url)

        with open(target + '.pdb.gz', 'wb') as fd:
            for chunk in r.iter_content(chunk_size=128):
                fd.write(chunk)

    if not os.path.isfile(target + '.pdb'):
        with gzip.open(target + '.pdb.gz', 'rb') as f_in:
            with open(target + '.pdb', 'wb') as f_out:
                f_out.write(f_in.read())
    else:
        print('[INFO] Pdb file %s.pdb already exists. Not downloading' % (
            target
        ))

    return target + '.pdb'


def load_3dab_model(filename):
    print('[INFO] Parsing %s' % (filename))

    data = []

    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            amino, posx, posy, posz = row
            data.append([
                amino,
                float(posx),
                float(posy),
                float(posz)]
            )

    return data


def main():
    csv_filename = sys.argv[1]
    target = sys.argv[2]

    filename = fetch_native_pdb(target)
    ab_model = load_3dab_model(csv_filename)

    from distance_normalizer import normalize_distances

    normalize_distances(filename, ab_model)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('usage: %s 3bAB_file.csv pdbcode' % sys.argv[0])
    else:
        main()

