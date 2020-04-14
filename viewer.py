#!/usr/bin/python3

import os, sys, curses, time, configparser, re, bisect, shutil

nucleotide_decoding = {
    0 : 'A',
    1 : 'C',
    2 : 'G',
    3 : 'T',
    4 : '?'
}

feature_encode = {
    'gap' : 0,
    'exon' : 1,
    'pseudogene' : 2,
    'gene' : 3
} # 0-63
feature_mask = 63
end_encode = 128

PAIR_UNK = 0
PAIR_NONE = 8
PAIR_EXON_PSEUDO = 12
PAIR_EXON_GENE = 16

nucleotide_colors = {
    0 : 9,
    1 : 11,
    2 : 10,
    3 : 14,
    4 : 5
}

region_colors = {
    PAIR_NONE : -1,
    PAIR_EXON_PSEUDO : 102,
    PAIR_EXON_GENE : 170,
}

script_path = os.path.realpath(__file__)
path = os.path.dirname(script_path)

scrx = 0
scry = 0
scrw = 80
scrh = 25

current_features = {}
current_ch = "1"
pos = 1
pos_percent = False
size = None
next_pos = None
next_feat = None

def get_input(stdscr):
    global scrw, scrh
    time.sleep(0.1)
    key = stdscr.getch()
    if key == curses.KEY_RESIZE:
        scrh, scrw = stdscr.getmaxyx()
        scrx, scry = 0, 0

def print_status(stdscr):
    global pos, size
    stdscr.addstr(0, 0, "{} ({:.3f}%)".format(pos, pos*100/size))

def next_line(stdscr):
    global scrx, scry, scrw, scrh
    scrx = 0
    scry += 1
    if scry >= scrh:
        stdscr.scroll()
        scry = scrh - 1
        print_status(stdscr)
        get_input(stdscr)

def next_char(stdscr):
    global scrx, scry, scrw, scrh, pos, size
    scrx += 1
    if scrx >= scrw-1:
        next_line(stdscr)

def print_char(char, pair, stdscr):
    global scrx, scry
    try:
        stdscr.addch(scry, scrx, char, curses.color_pair(pair))
        next_char(stdscr)
    except curses.error:
        next_char(stdscr)
        print_char(char, pair, stdscr)

def print_nucleotide(nucleotide, stdscr):
    global current_feature
    pair = None
    if feature_encode['gap'] in current_features:
        nucleotide = 4
        pair = PAIR_UNK
    elif feature_encode['exon'] in current_features:
        if feature_encode['gene'] in current_features:
            pair = PAIR_EXON_GENE + nucleotide
        elif feature_encode['pseudogene'] in current_features:
            pair = PAIR_EXON_PSEUDO + nucleotide
        else:
            pair = PAIR_UNK
    else:
        pair = PAIR_NONE + nucleotide
    print_char(nucleotide_decoding[nucleotide], pair, stdscr)

  #####      ##    #######   #######  ##        ########  #######  ########  #######   #######  ##     ## ##    ##
 ##   ##   ####   ##     ## ##     ## ##    ##  ##       ##     ## ##    ## ##     ## ##     ##  ##   ##   ##  ##                ##
##     ##    ##          ##        ## ##    ##  ##       ##            ##   ##     ## ##     ##   ## ##     ####   ## ##  ##  ########
##     ##    ##    #######   #######  ##    ##  #######  ########     ##     #######   ########    ###       ##    ### ### ##    ##
##     ##    ##   ##               ## #########       ## ##     ##   ##     ##     ##        ##   ## ##      ##    ##  ##  ##    ##
 ##   ##     ##   ##        ##     ##       ##  ##    ## ##     ##   ##     ##     ## ##     ##  ##   ##     ##    ##  ##  ##    ##
  #####    ###### #########  #######        ##   ######   #######    ##      #######   #######  ##     ##    ##    ##  ##  ##     ####

def print_title(title, stdscr):
    global scrw
    chars = {
        '0' : [
            "  #####    ",
            " ##   ##   ",
            "##     ##  ",
            "##     ##  ",
            "##     ##  ",
            " ##   ##   ",
            "  #####    ",
        ],
        '1' : [
            "  ##   ",
            "####   ",
            "  ##   ",
            "  ##   ",
            "  ##   ",
            "  ##   ",
            "###### ",
        ],
        '2' : [
            " #######  ",
            "##     ## ",
            "       ## ",
            " #######  ",
            "##        ",
            "##        ",
            "######### ",
        ],
        '3' : [
            " #######  ",
            "##     ## ",
            "       ## ",
            " #######  ",
            "       ## ",
            "##     ## ",
            " #######  ",
        ],
        '4' : [
            "##        ",
            "##    ##  ",
            "##    ##  ",
            "##    ##  ",
            "######### ",
            "      ##  ",
            "      ##  ",
        ],
        '5' : [
            "######## ",
            "##       ",
            "##       ",
            "#######  ",
            "      ## ",
            "##    ## ",
            " ######  ",
        ],
        '6' : [
            " #######  ",
            "##     ## ",
            "##        ",
            "########  ",
            "##     ## ",
            "##     ## ",
            " #######  ",
        ],
        '7' : [
            "######## ",
            "##    ## ",
            "    ##   ",
            "   ##    ",
            "  ##     ",
            "  ##     ",
            "  ##     ",
        ],
        '8' : [
            " #######  ",
            "##     ## ",
            "##     ## ",
            " #######  ",
            "##     ## ",
            "##     ## ",
            " #######  ",
        ],
        '9' : [
            " #######  ",
            "##     ## ",
            "##     ## ",
            " ######## ",
            "       ## ",
            "##     ## ",
            " #######  ",
        ],
        'X' : [
            "##     ## ",
            " ##   ##  ",
            "  ## ##   ",
            "   ###    ",
            "  ## ##   ",
            " ##   ##  ",
            "##     ## ",
        ],
        'Y' : [
            "##    ## ",
            " ##  ##  ",
            "  ####   ",
            "   ##    ",
            "   ##    ",
            "   ##    ",
            "   ##    ",
        ],
        'm' : [
            "           ",
            "           ",
            "## ##  ##  ",
            "### ### ## ",
            "##  ##  ## ",
            "##  ##  ## ",
            "##  ##  ## ",
        ],
        't' : [
            "         ",
            "   ##    ",
            "######## ",
            "   ##    ",
            "   ##    ",
            "   ##    ",
            "    #### ",
        ],
    }
    length = 0
    for C in title:
        length += len(chars[C][0])
    next_line(stdscr)
    pad = (scrw - length)//2
    if pad < 0:
        pad = 0
    for n in range(0, len(chars['0'])):
        for m in range(0, pad):
            print_char(" ", PAIR_UNK, stdscr)
        for C in title:
            for c in chars[C][n]:
                print_char(c, PAIR_UNK, stdscr)
        next_line(stdscr)
    next_line(stdscr)

def update_features(mt_file):
    global next_pos, next_feat, current_features
    if next_feat & end_encode:
        if (next_feat & feature_mask) not in current_features:
            pass
        current_features[next_feat & feature_mask] -= 1
        if current_features[next_feat & feature_mask] == 0:
            del current_features[next_feat & feature_mask]
    else:
        if (next_feat & feature_mask) not in current_features:
            current_features[next_feat & feature_mask] = 0
        current_features[next_feat & feature_mask] += 1
    next_pos = int.from_bytes(mt_file.read(4), byteorder='little', signed=False)
    next_feat = int.from_bytes(mt_file.read(1), byteorder='little', signed=False)

def get_start_pos():
    global current_ch, pos, pos_percent
    match = None
    for arg in sys.argv[1:]:
        match = re.fullmatch(r'([1-9XY]|1\d|2[0-2]|mt)\s*(\.(-?\d+%?))?', arg)
        if match:
            break
    if match:
        current_ch = match.group(1)
        pos_str = match.group(3)
        if pos_str:
            if pos_str[-1] == '%':
                pos_percent = True
                pos = int(pos_str[:-1])
            else:
                pos = int(pos_str)

def main(stdscr):
    global next_pos, next_feat, current_features, pos, size, pos_percent
    curses.start_color()
    curses.use_default_colors()
    stdscr.idlok(True)
    stdscr.scrollok(True)
    stdscr.immedok(False)
    stdscr.nodelay(True)
    for pair, background in region_colors.items():
        for offset, foreground in nucleotide_colors.items():
            curses.init_pair(pair + offset, foreground, background)

    ch_path = os.path.join(path, current_ch + ".bin")
    file = open(ch_path, 'rb')
    size = file.read(4)
    size = int.from_bytes(size, byteorder='little', signed=False)

    if pos_percent:
        pos = (pos * size) // 100
    if pos <= 0:
        pos = size + pos
    file.seek(4 + (pos-1)//4)

    mt_path = os.path.join(path, current_ch + ".dat")
    mt_file = open(mt_path, 'rb')
    next_pos = int.from_bytes(mt_file.read(4), byteorder='little', signed=False)
    next_feat = int.from_bytes(mt_file.read(1), byteorder='little', signed=False)

    while(pos > next_pos and next_pos > 0):
        update_features(mt_file)

    byte = file.read(1)
    while byte != b"":
        byte = int.from_bytes(byte, byteorder='little')
        for n in range(0, 4):
            if pos == 1:
                for N in range(0, 10):
                    next_line(stdscr)
                print_title(current_ch, stdscr)
            if pos == size:
                #get_input(stdscr)
                #sleep(10.0)
                break
            if pos == next_pos:
                while pos == next_pos:
                    update_features(mt_file)
            pos += 1
            nucleotide = (byte >> (2*(3-n))) & 3
            print_nucleotide(nucleotide, stdscr)
        byte = file.read(1)
    stdscr.getch()

def index_closest(list, val):
    pos = bisect.bisect_left(list, val)
    if pos in (0, len(list)):
        return pos
    before = list[pos - 1]
    after = list[pos]
    if after - val < val - before:
       return pos
    else:
       return pos - 1

def convert_color(R, G, B):
    basic = {
        (000, 000, 000) : 0,
        (128, 000, 000) : 1,
        (000, 120, 000) : 2,
        (128, 128, 000) : 3,
        (000, 000, 128) : 4,
        (128, 000, 128) : 5,
        (000, 128, 128) : 6,
        (192, 192, 192) : 7,
        (128, 128, 128) : 8
        #the others can be mapped later
    }
    if (R, G, B) in basic:
        return basic[(R, G, B)]

    if R == G == B and R < 243:
        V = (R - 3)//10
        if V < 0:
            V = 0
        return V + 232

    steps = [0, 95, 135, 175, 215, 255]
    R = index_closest(steps, R)
    G = index_closest(steps, G)
    B = index_closest(steps, B)
    return R*36 + G*6 + B + 16

def get_config_color(dict, key, section, option):
    str = section.get(option)
    if not str:
        return
    if re.fullmatch(r'\d+', str) and int(str) in range(0, 256):
        dict[key] = int(str)
        return
    match = re.fullmatch(r'#([\da-fA-F]{2})([\da-fA-F]{2})([\da-fA-F]{2})', str)
    if match:
        R = int(match.group(1), base=16)
        G = int(match.group(2), base=16)
        B = int(match.group(3), base=16)
        dict[key] = convert_color(R, G, B)
        return
    match = re.fullmatch(r'[rR][gG][bB]\((\d+),\s*(\d+),\s*(\d+)\)', str)
    if match:
        R = int(match.group(1), base=10)
        G = int(match.group(2), base=10)
        B = int(match.group(3), base=10)
        dict[key] = convert_color(R, G, B)
        return

def parse_config():
    config = configparser.ConfigParser()
    config_path = os.path.join(path, "config.ini")
    config.read(config_path)
    if 'Nucleobase Colors' in config:
        section = config['Nucleobase Colors']
        get_config_color(nucleotide_colors, 0, section, 'A')
        get_config_color(nucleotide_colors, 1, section, 'C')
        get_config_color(nucleotide_colors, 2, section, 'G')
        get_config_color(nucleotide_colors, 3, section, 'T')
        get_config_color(nucleotide_colors, 4, section, '?')
    if 'Region Colors' in config:
        section = config['Region Colors']
        get_config_color(region_colors, PAIR_EXON_PSEUDO, section, 'pseudogene exon')
        get_config_color(region_colors, PAIR_EXON_GENE, section, 'gene exon')

parse_config()
get_start_pos()
scrw, scrh = shutil.get_terminal_size((scrw, scrh))
curses.wrapper(main)
