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
    'CDS' : 2,
    'pseudogene' : 3,
    'gene' : 4,
    'tRNA' : 5,
    'rRNA' : 6,
    'miRNA' : 7
} # 0-63
feature_mask = 63
end_encode = 128

strand_decode = {
    0 : '.',
    1 : '+',
    2 : '-'
}

PAIR_UNK = 0
PAIR_HIGHLIGHT = 1
PAIR_NONE = 8
PAIR_EXON_PSEUDO = 12
PAIR_UTR_GENE = 16
PAIR_CDS = 20
PAIR_CDS2 = 24
PAIR_INTRON = 28
PAIR_TRNA = 32
PAIR_RRNA = 36
PAIR_MIRNA = 40

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
    PAIR_UTR_GENE : 170,
    PAIR_CDS : 63,
    PAIR_CDS2 : 105,
    PAIR_INTRON : 232,
    PAIR_TRNA : 106,
    PAIR_RRNA : 65,
    PAIR_MIRNA : 136
}

other_colors = {
    PAIR_HIGHLIGHT : 11
}

script_path = os.path.realpath(__file__)
path = os.path.dirname(script_path)

scrx = 0
scry = 0
scrw = 80
scrh = 25

highlight = {
    'cpg' : False
}

current_info_strand = 0
ch_initial = "1"
pos_initial = 1
pos_percent = False
prev_nucleotide = None
paused = False
current_frame = 0

class Reader:
    def get_feature_info(self):
        global current_info_strand, current_frame
        if self.next_feat == feature_encode['gene']:
            info = b""
            current_info_strand = int.from_bytes(self.mt_file.read(1), byteorder='little')
            byte = self.mt_file.read(1)
            while byte != b"\0":
                info += byte
                byte = self.mt_file.read(1)
            return info.decode()
        if self.next_feat == feature_encode['CDS']:
            phase = int.from_bytes(self.mt_file.read(1), byteorder='little')
            current_frame = (current_frame & 4) | (3 - phase) % 3
        return ""

    def update_features(self):
        if self.next_feat & end_encode:
            if (self.next_feat & feature_mask) not in self.current_features:
                pass
            self.current_features[self.next_feat & feature_mask] -= 1
            if self.current_features[self.next_feat & feature_mask] == 0:
                del self.current_features[self.next_feat & feature_mask]
        else:
            if (self.next_feat & feature_mask) not in self.current_features:
                self.current_features[self.next_feat & feature_mask] = 0
            self.current_features[self.next_feat & feature_mask] += 1
        info = self.get_feature_info()
        if info:
            self.current_info = info
            self.prev_info_pos = self.next_pos
        self.next_pos = int.from_bytes(self.mt_file.read(4), byteorder='little', signed=False)
        self.next_feat = int.from_bytes(self.mt_file.read(1), byteorder='little', signed=False)

    def __init__(self):
        global ch_initial, pos_initial, pos_percent

        self.ch = ch_initial
        ch_path = os.path.join(path, self.ch + ".bin")
        self.file = open(ch_path, 'rb')
        ch_size = self.file.read(4)
        self.ch_size = int.from_bytes(ch_size, byteorder='little', signed=False)

        self.pos = pos_initial
        if pos_percent:
            self.pos = (self.pos * self.ch_size) // 100
        if pos_initial <= 0:
            self.pos = self.ch_size + self.pos
        self.file.seek(4 + (self.pos-1)//4)

        mt_path = os.path.join(path, self.ch + ".dat")
        self.mt_file = open(mt_path, 'rb')
        self.next_pos = int.from_bytes(self.mt_file.read(4), byteorder='little', signed=False)
        self.next_feat = int.from_bytes(self.mt_file.read(1), byteorder='little', signed=False)

        self.current_features = {}
        self.current_info = ""
        self.prev_info_pos = None

        while(self.pos > self.next_pos and self.next_pos > 0):
            self.update_features()

        self.byte = self.file.read(1)
        self.n = (self.pos-1) % 4
        if self.byte == b"":
            self.eof = True
        else:
            self.eof = False

    def read(self):
        b = int.from_bytes(self.byte, byteorder='little')
        return (b >> (2*(3-self.n))) & 3

    def advance(self):
        global scrw, scrh
        self.n += 1
        if self.n & 3 == 0:
            self.byte = self.file.read(1)
            self.n &= 3
        self.pos += 1
        if self.byte == b"":
            self.eof = True
        if self.pos == self.next_pos:
            while self.pos == self.next_pos:
                self.update_features()
        if self.current_info and self.prev_info_pos and self.pos - self.prev_info_pos > scrw * scrh:
            self.current_info = ""

    def __del__(self):
        self.file.close()
        self.mt_file.close()

def get_input(stdscr):
    global scrw, scrh, paused
    time.sleep(0.1)
    key = stdscr.getch()
    if key == curses.KEY_RESIZE:
        scrh, scrw = stdscr.getmaxyx()
        scrx, scry = 0, 0
    elif key == ord('\n') or key == curses.KEY_ENTER or key == ord(' '):
        paused = not paused

class View:
    def __init__(self, reader, stdscr):
        self.reader = reader
        self.screen = stdscr

    def print_status(self):
        status = "{} ({:.3f}%)".format(self.reader.pos, self.reader.pos*100/self.reader.ch_size)
        if self.reader.current_info:
            status += " {} ({})".format(self.reader.current_info, strand_decode[current_info_strand])

        self.screen.addstr(0, 0, status)

    def next_line(self):
        global scrx, scry, scrw, scrh, paused
        scrx = 0
        scry += 1
        if scry >= scrh:
            self.print_status()
            get_input(self.screen)
            while(paused):
                get_input(self.screen)
            self.screen.scroll()
            scry = scrh - 1

    def next_char(self):
        global scrx, scry, scrw, scrh, pos, size
        scrx += 1
        if scrx >= scrw-1:
            self.next_line()

    def print_char(self, char, pair):
        global scrx, scry
        try:
            self.screen.addch(scry, scrx, char, curses.color_pair(pair))
            self.next_char()
        except curses.error:
            self.next_char()
            self.print_char(char, pair)

    def set_prev_pairs(self, number, pair):
        global scrx, scry, scrw, scrh
        x = scrx
        y = scry
        for n in range(0, number):
            x -= 1
            if x < 0:
                x = scrw-2
                y -= 1
                if y < 0:
                    break
            self.screen.chgat(y, x, curses.color_pair(pair))

    def print_nucleotide(self):
        global prev_nucleotide, highlight, current_frame
        pair = None
        nucleotide = self.reader.read()
        features = self.reader.current_features
        if feature_encode['gap'] in features:
            nucleotide = 4
            pair = PAIR_UNK
        elif feature_encode['CDS'] in features:
            if current_frame & 4:
                pair = PAIR_CDS2 + nucleotide
            else:
                pair = PAIR_CDS + nucleotide
            current_frame += 1
            if current_frame & 3 == 3:
                current_frame ^= 4
                current_frame &= 4
        elif feature_encode['tRNA'] in features:
            pair = PAIR_TRNA + nucleotide
        elif feature_encode['rRNA'] in features:
            pair = PAIR_RRNA + nucleotide
        elif feature_encode['miRNA'] in features:
            pair = PAIR_MIRNA + nucleotide
        elif feature_encode['exon'] in features:
            if feature_encode['gene'] in features:
                pair = PAIR_UTR_GENE + nucleotide
            elif feature_encode['pseudogene'] in features:
                pair = PAIR_EXON_PSEUDO + nucleotide
            else:
                pair = PAIR_UNK
        elif feature_encode['gene'] in features or feature_encode['pseudogene'] in features:
            pair = PAIR_INTRON + nucleotide
        else:
            pair = PAIR_NONE + nucleotide
        if highlight['cpg'] and prev_nucleotide == 1 and nucleotide == 2:
            pair = PAIR_HIGHLIGHT
            self.set_prev_pairs(1, PAIR_HIGHLIGHT)
        prev_nucleotide = nucleotide
        self.print_char(nucleotide_decoding[nucleotide], pair)

    def print_title(self, title):
          #####      ##    #######   #######  ##        ########  #######  ########  #######   #######  ##     ## ##    ##
         ##   ##   ####   ##     ## ##     ## ##    ##  ##       ##     ## ##    ## ##     ## ##     ##  ##   ##   ##  ##                ##
        ##     ##    ##          ##        ## ##    ##  ##       ##            ##   ##     ## ##     ##   ## ##     ####   ## ##  ##  ########
        ##     ##    ##    #######   #######  ##    ##  #######  ########     ##     #######   ########    ###       ##    ### ### ##    ##
        ##     ##    ##   ##               ## #########       ## ##     ##   ##     ##     ##        ##   ## ##      ##    ##  ##  ##    ##
         ##   ##     ##   ##        ##     ##       ##  ##    ## ##     ##   ##     ##     ## ##     ##  ##   ##     ##    ##  ##  ##    ##
          #####    ###### #########  #######        ##   ######   #######    ##      #######   #######  ##     ##    ##    ##  ##  ##     ####

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
        self.next_line()
        pad = (scrw - length)//2
        if pad < 0:
            pad = 0
        for n in range(0, len(chars['0'])):
            for m in range(0, pad):
                self.print_char(" ", PAIR_UNK)
            for C in title:
                for c in chars[C][n]:
                    self.print_char(c, PAIR_UNK)
            self.next_line()
        self.next_line()

    def scroll_down(self, n):
        while n > 0 and not self.reader.eof:
            if self.reader.pos == 1:
                for N in range(0, 10):
                    self.next_line()
                self.print_title(self.reader.ch)
            if self.reader.pos == self.reader.ch_size:
                break
            self.print_nucleotide()
            self.reader.advance()

    def __del__(self):
        pass

def main(stdscr):
    global paused, current_reader
    curses.start_color()
    curses.use_default_colors()
    stdscr.idlok(True)
    stdscr.scrollok(True)
    stdscr.immedok(False)
    stdscr.nodelay(True)
    for pair, background in region_colors.items():
        for offset, foreground in nucleotide_colors.items():
            curses.init_pair(pair + offset, foreground, background)
    curses.init_pair(PAIR_HIGHLIGHT, 0, other_colors[PAIR_HIGHLIGHT])

    reader = Reader()
    view = View(reader, stdscr)

    while not reader.eof:
        view.scroll_down(1)

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
        get_config_color(region_colors, PAIR_UTR_GENE, section, 'gene UTR')
        get_config_color(region_colors, PAIR_CDS, section, 'CDS')
        get_config_color(region_colors, PAIR_CDS2, section, 'CDS 2')
        get_config_color(region_colors, PAIR_INTRON, section, 'intron')
        get_config_color(region_colors, PAIR_TRNA, section, 'tRNA')
        get_config_color(region_colors, PAIR_RRNA, section, 'rRNA')
        get_config_color(region_colors, PAIR_MIRNA, section, 'miRNA')
    if 'Other Colors' in config:
        get_config_color(other_colors, PAIR_HIGHLIGHT, section, 'highlight')

def get_start_pos():
    global ch_initial, pos_initial, pos_percent, paused
    match = None
    for arg in sys.argv[1:]:
        match = re.fullmatch(r'([1-9XY]|1\d|2[0-2]|mt)\s*(\.(-?\d+%?))?', arg)
        if match:
            break
    if match:
        ch_initial = match.group(1)
        pos_str = match.group(3)
        if pos_str:
            if pos_str[-1] == '%':
                pos_percent = True
                pos_initial = int(pos_str[:-1])
            else:
                pos_initial = int(pos_str)
            paused = True

def parse_options():
    global highlight
    get_start_pos()
    match = None
    for arg in sys.argv[1:]:
        match = re.fullmatch(r'hl=([a-zA-Z0-9,]*)', arg)
        if match:
            break
    if match:
        hls = match.group(1).split(',')
        for hl in hls:
            hl = hl.lower()
            if hl in highlight:
                highlight[hl] = True

parse_config()
parse_options()
scrw, scrh = shutil.get_terminal_size((scrw, scrh))
curses.wrapper(main)
