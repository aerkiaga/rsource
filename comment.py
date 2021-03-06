#!/usr/bin/pypy3

import sys, os, re, bisect

feature_encode = {
    'gap' : 0,
    'exon' : 1,
    'CDS' : 2,
    'pseudogene' : 3,
    'gene' : 4,
    'tRNA' : 5,
    'rRNA' : 6,
    'miRNA' : 7
} # 0-63
end_encode = 128

strand_encode = {
    '+' : 1,
    '-' : 2,
    '.' : 3
}

script_path = os.path.realpath(__file__)
path = os.path.dirname(script_path)

current_ch = None
ch_files = {}
ch_arr_pos = {}
ch_arr_feat = {}
ch_arr_info = {}

pattern_seqid = re.compile(r'NC_(\d+)')
pattern_info_description = re.compile(r';description=([^;]*);')
pattern_info_name = re.compile(r';Name=([^;]*);')

def insert_feature(pos, feat, info=None):
    index = bisect.bisect_right(ch_arr_pos[current_ch], pos)
    ch_arr_pos[current_ch].insert(index, pos)
    ch_arr_feat[current_ch].insert(index, feat)
    ch_arr_info[current_ch].insert(index, info)

def get_feature_info(feat, fields):
    if feat == feature_encode['gene']:
        match = pattern_info_description.search(fields[8])
        if not match:
            match = pattern_info_name.search(fields[8])
        if not match:
            return None
        info = b'\0'
        info += strand_encode[fields[6]].to_bytes(1, byteorder='little')
        info += match.group(1).encode()
        info += b'\0'
        return info
    elif feat == feature_encode['CDS']:
        info = int(fields[7]).to_bytes(1, byteorder='little')
        return info
    return None

for line in sys.stdin:
    if line[0] != '#':
        fields = line.split('\t')
        match = pattern_seqid.match(fields[0])
        if not match:
            continue
        if fields[2] not in feature_encode:
            continue

        current_ch = int(match.group(1))
        if current_ch == 23:
            current_ch = 'X'
        elif current_ch == 24:
            current_ch = 'Y'
        elif current_ch == 12920:
            current_ch = 'mt'
        else:
            current_ch = str(current_ch)

        if current_ch not in ch_files:
            current_ch_path = os.path.join(path, current_ch + ".dat")
            ch_files[current_ch] = open(current_ch_path, 'wb')
            ch_arr_pos[current_ch] = []
            ch_arr_feat[current_ch] = []
            ch_arr_info[current_ch] = []
            if current_ch == 'mt':
                print("Mitochondrial")
            else:
                print("Chromosome " + current_ch)

        pos = int(fields[3])
        endpos = int(fields[4])
        endpos += 1
        feat = feature_encode[fields[2]]
        info = get_feature_info(feat, fields)
        insert_feature(pos, feat, info)
        insert_feature(endpos, feat | end_encode)

print("Annotating gaps and saving...")
for ch in ch_files.keys():
    current_ch = ch
    gap_file_path = os.path.join(path, ch + ".gap")
    gap_file = open(gap_file_path, 'rb')
    gap_start = gap_file.read(4)
    gap_end = gap_file.read(4)
    while gap_start != b"":
        gap_start = int.from_bytes(gap_start, byteorder='little', signed=False)
        gap_end = int.from_bytes(gap_end, byteorder='little', signed=False)
        if gap_end != 0:
            insert_feature(gap_start, feature_encode['gap'])
            insert_feature(gap_end, feature_encode['gap'] | end_encode)
        gap_start = gap_file.read(4)
        gap_end = gap_file.read(4)

    file = ch_files[ch]
    for n in range(0, len(ch_arr_pos[ch])):
        file.write(ch_arr_pos[ch][n].to_bytes(4, byteorder='little', signed=False))
        file.write(ch_arr_feat[ch][n].to_bytes(1, byteorder='little', signed=False))
        if ch_arr_info[ch][n]:
            file.write(ch_arr_info[ch][n])
            file.write(ch_arr_feat[ch][n].to_bytes(1, byteorder='little', signed=False))
print("Done!")
