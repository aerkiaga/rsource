#!/usr/bin/pypy3

import sys, os, re

chromosome_strings = {
    '1' : "249 Mbp ",
    '2' : "242 Mbp      ",
    '3' : "198 Mbp      ",
    '4' : "190 Mbp             ",
    '5' : "182 Mbp             ",
    '6' : "171 Mbp          ",
    '7' : "159 Mbp           ",
    '8' : "145 Mbp             ",
    '9' : "138 Mbp             ",
    '10' : "134 Mbp             ",
    '11' : "135 Mbp           ",
    '12' : "133 Mbp              ",
    '13' : "114 Mbp                ",
    '14' : "107 Mbp                ",
    '15' : "102 Mbp                ",
    '16' : " 90 Mbp              ",
    '17' : " 83 Mbp               ",
    '18' : " 80 Mbp                 ",
    '19' : " 59 Mbp               ",
    '20' : " 64 Mbp               ",
    '21' : " 47 Mbp                ",
    '22' : " 51 Mbp                ",
    'X' : "156 Mbp          ",
    'Y' : " 57 Mbp                  ",
    'mt' : " 17 kbp                    "
}

chromosome_lengths = {
    '1' : 248956422,
    '2' : 242193529,
    '3' : 198295559,
    '4' : 190214555,
    '5' : 181538259,
    '6' : 170805979,
    '7' : 159345973,
    '8' : 145138636,
    '9' : 138394717,
    '10' : 133797422,
    '11' : 135086622,
    '12' : 133275309,
    '13' : 114364328,
    '14' : 107043718,
    '15' : 101991189,
    '16' : 90338345,
    '17' : 83257441,
    '18' : 80373285,
    '19' : 58617616,
    '20' : 64444167,
    '21' : 46709983,
    '22' : 50818468,
    'X' : 156040895,
    'Y' : 57227415,
    'mt' : 16569
}

chromosome_progress = {
    '1' : "(█▐█ ▌▌▐▌___▐█▌__██╳▒▒▒██_█_█___██▌_█▌▐)",
    '2' : "(█_█▌_█▌_▌█▌_█╳█_▌▐▐█__█__█▌__█▌▐▌▐██)",
    '3' : "(▐▌_▐_███__█__╳▒_▌_█▌▌_█_▐_█_█)",
    '4' : "(█▌▐_█▐╳█_▐█_▌▐_▌__█▐▌▐_█)",
    '5' : "(█▐__█▌╳█_██__▌_▐__██▌▐▌_█)",
    '6' : "(█▐▌_█▌█_▐╳_▌▐▌_██___█__█▐█)",
    '7' : "(█__██_█▐╳██▌__██___█▌▌██)",
    '8' : "(█▌▐█_▐╳▌▐█_▐_▐█▌_▐█▐█)",
    '9' : "(█_▌_█_╳▒▒▌__██▌_▌▐██)",
    '10' : "(█_█_▐█╳█___██▌_█▌_██)",
    '11' : "(███__█_█╳_██▌_▐_▐██▌█)",
    '12' : "(██__█╳_██_▐▌__█_█▐█)",
    '13' : "(▒═▒╳█▌▐█▌__█▌_█_█)",
    '14' : "(▒═▒╳█_█__█▌██▌_█)",
    '15' : "(▒═▒╳▌▌▐▌_▐█▐█_██)",
    '16' : "(█▐▌_█╳▒▌█_█▌▐█)",
    '17' : "(█▌_█╳█_██_▐_▐█)",
    '18' : "(_█╳█__██__█)",
    '19' : "(█▌▐▒╳▒█▌▐█_)",
    '20' : "(█_▐█╳█▌▐▌█)",
    '21' : "(▒═▒╳▌_██)",
    '22' : "(▒═▒╳█▌▐██)",
    'X' : "(█▌█▌_▐▌██╳▐██___█▌▐▌_█▌▐█)",
    'Y' : "(▐╳██▒▒▒)",
    'mt' : "o"
}

nucleotide_encoding = {
    'A' : 0,
    'C' : 1,
    'G' : 2,
    'T' : 3,
}

wildcard_nucleotides = ['B', 'D', 'H', 'K', 'M', 'N', 'R', 'S', 'V', 'W', 'Y']

script_path = os.path.realpath(__file__)
path = os.path.dirname(script_path)

current_ch = None
last_ch = None
ch_files = {}
ch_lengths = {}
ch_bytes = {}
ch_progress = {}

gap_files = {}
gap_starts = {}

pattern_chromosome = re.compile(r'>.+?Homo sapiens chromosome ([1-9XY]|1\d|2[0-2]), GRCh.+?Primary Assembly')
pattern_mitichondrial = re.compile(r'>.+?Homo sapiens mitochondrion, complete genome')

def close_current_chromosome(ch):
    if not ch:
        return
    if ch_lengths[ch] % 4:
        ch_bytes[ch] <<= (8 - 2*(ch_lengths[ch] % 4))
        ch_files[ch].write(ch_bytes[ch].to_bytes(1, byteorder='little', signed=False))
    ch_files[ch].seek(0)
    ch_files[ch].write(ch_lengths[ch].to_bytes(4, byteorder='little', signed=False))
    ch_files[ch].close()
    if ch in gap_starts:
        gap_files[ch].write(ch_lengths[ch].to_bytes(4, byteorder='little', signed=False))
        del gap_starts[ch]
    gap_files[ch].close()
    if ch_progress[ch] < len(chromosome_progress[ch]) - 1:
        print(chromosome_progress[ch][ch_progress[ch]+1:], end='', flush=True)

for line in sys.stdin:
    if(line[0] == '>' or line[0] == ';'):
        match = pattern_chromosome.match(line)
        if match:
            close_current_chromosome(last_ch)
            last_ch = current_ch = match.group(1)
            print("\nChromosome " + current_ch + "\t" + chromosome_strings[current_ch], end='')
        else:
            match = pattern_mitichondrial.match(line)
            if match:
                close_current_chromosome(last_ch)
                last_ch = current_ch = 'mt'
                print("\nMitochondrial\t" + chromosome_strings[current_ch], end='')
            else:
                current_ch = None
    else:
        if current_ch:
            if current_ch not in ch_files:
                current_ch_path = os.path.join(path, current_ch + ".bin")
                ch_files[current_ch] = open(current_ch_path, 'wb')
                ch_files[current_ch].write((0).to_bytes(4, byteorder='little', signed=False))
                ch_lengths[current_ch] = 0
                ch_bytes[current_ch] = 0
                ch_progress[current_ch] = -1
            ba = bytearray()
            line = line.upper()
            for c in line:
                if c in nucleotide_encoding:
                    ch_bytes[current_ch] = (ch_bytes[current_ch] << 2) | nucleotide_encoding[c]
                    ch_lengths[current_ch] += 1
                    if ch_lengths[current_ch] % 4 == 0:
                        ba.append(ch_bytes[current_ch])
                        ch_bytes[current_ch] = 0
                    if current_ch in gap_starts:
                        gap_files[current_ch].write(ch_lengths[current_ch].to_bytes(4, byteorder='little', signed=False))
                        del gap_starts[current_ch]
                elif c in wildcard_nucleotides:
                    ch_bytes[current_ch] = (ch_bytes[current_ch] << 2)
                    ch_lengths[current_ch] += 1
                    if ch_lengths[current_ch] % 4 == 0:
                        ba.append(ch_bytes[current_ch])
                        ch_bytes[current_ch] = 0
                    if current_ch not in gap_starts:
                        gap_starts[current_ch] = ch_lengths[current_ch]
                        if current_ch not in gap_files:
                            current_gap_path = os.path.join(path, current_ch + ".gap")
                            gap_files[current_ch] = open(current_gap_path, 'wb')
                        gap_files[current_ch].write(gap_starts[current_ch].to_bytes(4, byteorder='little', signed=False))

            ch_files[current_ch].write(ba)
            progress = ch_lengths[current_ch] * len(chromosome_progress[current_ch]) // chromosome_lengths[current_ch]
            if(progress > ch_progress[current_ch]):
                if progress < len(chromosome_progress[current_ch]):
                    print(chromosome_progress[current_ch][progress], end='', flush=True)
                ch_progress[current_ch] += 1
close_current_chromosome(last_ch)
print("\nDone!")
