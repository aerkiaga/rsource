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
    1 : '+',
    2 : '-',
    3 : '.'
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

ch_initial = "1"
pos_initial = 1
pos_percent = False
prev_nucleotide = None
paused = False

class Reader:
    def apply_feature(self, feat):
        if feat & end_encode:
            if (feat & feature_mask) not in self.current_features:
                return #D
            self.current_features[feat & feature_mask] -= 1
            if self.current_features[feat & feature_mask] == 0:
                del self.current_features[feat & feature_mask]
        else:
            if (feat & feature_mask) not in self.current_features:
                self.current_features[feat & feature_mask] = 0
            self.current_features[feat & feature_mask] += 1

    def get_feature_info(self):
        if self.next_feat == feature_encode['gene']:
            info = b""
            self.mt_file.read(1)
            strand = self.mt_file.read(1)
            info += strand
            byte = self.mt_file.read(1)
            while byte != b"\0":
                info += byte
                byte = self.mt_file.read(1)
            self.mt_file.read(1)
            return info.decode()
        if self.next_feat == feature_encode['CDS']:
            phase = int.from_bytes(self.mt_file.read(1), byteorder='little')
            self.mt_file.read(1)
        return ""

    def update_features(self):
        self.cur_feat_pos = self.next_pos
        self.apply_feature(self.next_feat)
        info = self.get_feature_info()
        if info:
            self.current_info_strand = ord(info[0])
            self.current_info = info[1:]
            self.prev_info_pos = self.next_pos
        self.next_pos = int.from_bytes(self.mt_file.read(4), byteorder='little', signed=False)
        self.next_feat = int.from_bytes(self.mt_file.read(1), byteorder='little', signed=False)

    def unget_feature(self):
        #only seeks, no other side effect; returns type of feature it lands in
        if self.mt_file.tell() < 10:
            return None
        self.mt_file.seek(self.mt_file.tell() - 6)
        cur_feat = int.from_bytes(self.mt_file.read(1), byteorder='little', signed=False)
        if cur_feat == feature_encode['gene']:
            self.mt_file.seek(self.mt_file.tell() - 3)
            while self.mt_file.read(1) != b'\0':
                self.mt_file.seek(self.mt_file.tell() - 2)
            self.mt_file.seek(self.mt_file.tell() - 1)
        elif cur_feat == feature_encode['CDS']:
            self.mt_file.seek(self.mt_file.tell() - 2)
        return cur_feat

    def update_features_backwards(self):
        self.next_pos = self.cur_feat_pos
        # we are after the new next-to-next's type
        lost_feat = self.unget_feature() # new next (previous current)
        self.apply_feature(lost_feat ^ end_encode)
        # we are after the new next's type
        exists = self.unget_feature() # new current
        if exists is not None:
            # we are after the new current's type
            self.mt_file.seek(self.mt_file.tell() - 5)
            # we are at the new current
            self.cur_feat_pos = int.from_bytes(self.mt_file.read(4), byteorder='little', signed=False) # position
            self.next_feat = int.from_bytes(self.mt_file.read(1), byteorder='little', signed=False) # type (actually of the current one)
            self.get_feature_info() # skip info with current type
            # we are at the next
            self.mt_file.read(5)
            # we are after the new next's type
        else:
            # we are after the new next's type
            self.cur_feat_pos = None
        self.next_feat = lost_feat # new next (previous current)

    def get_cds_phase(self):
        saved_fpos = self.mt_file.tell()
        if saved_fpos == self.cds_phase_cache['saved_fpos']:
            start_phase = self.cds_phase_cache['start_phase']
            start_pos = self.cds_phase_cache['start_pos']
            relative_pos = self.pos - start_pos
            r = (start_phase + relative_pos%3) | (((relative_pos // 3) & 1) << 2)
        else:
            prev_feat = self.unget_feature()
            while prev_feat != feature_encode['CDS'] and prev_feat is not None:
                prev_feat = self.unget_feature()
            if prev_feat == feature_encode['CDS']:
                self.cds_phase_cache['saved_fpos'] = saved_fpos
                self.cds_phase_cache['start_phase'] = start_phase = (3 - int.from_bytes(self.mt_file.read(1), byteorder='little', signed=False)) % 3
                self.mt_file.seek(self.mt_file.tell() - 6)
                self.cds_phase_cache['start_pos'] = start_pos = int.from_bytes(self.mt_file.read(4), byteorder='little', signed=False)
                relative_pos = self.pos - start_pos
                r = (start_phase + relative_pos%3) | (((relative_pos // 3) & 1) << 2)
            else:
                r = 0
            self.mt_file.seek(saved_fpos)
        return r

    def get_byte(self):
        self.byte = self.file.read(1)
        if self.byte == b'':
            self.eof = True
        else:
            self.eof = False

    def seek_pos(self):
        self.file.seek(4 + (self.pos-1)//4)
        self.n = (self.pos-1) % 4
        self.get_byte()

    def __init__(self, ch, pos, pos_is_percent):
        self.ch = ch
        ch_path = os.path.join(path, self.ch + ".bin")
        self.file = open(ch_path, 'rb')
        ch_size = self.file.read(4)
        self.ch_size = int.from_bytes(ch_size, byteorder='little', signed=False)

        self.pos = pos
        if pos_is_percent:
            self.pos = (self.pos * self.ch_size) // 100
        if pos_initial <= 0:
            self.pos = self.ch_size + self.pos
        self.seek_pos()

        mt_path = os.path.join(path, self.ch + ".dat")
        self.mt_file = open(mt_path, 'rb')
        self.next_pos = int.from_bytes(self.mt_file.read(4), byteorder='little', signed=False)
        self.next_feat = int.from_bytes(self.mt_file.read(1), byteorder='little', signed=False)

        self.current_features = {}
        self.cur_feat_pos = None
        self.current_info = ""
        self.prev_info_pos = None

        while(self.pos >= self.next_pos and self.next_pos > 0):
            self.update_features()

        self.cds_phase_cache = {
            'saved_fpos' : None,
            'start_phase' : None,
            'start_pos' : None
        }

    def read(self):
        b = int.from_bytes(self.byte, byteorder='little')
        return (b >> (2*(3-self.n))) & 3

    def advance(self):
        global scrw, scrh
        self.n += 1
        if self.n & 3 == 0:
            self.get_byte()
            self.n &= 3
        self.pos += 1
        if self.pos == self.next_pos:
            while self.pos == self.next_pos:
                self.update_features()

    def jump_to(self, P):
        if P > self.pos:
            while(P >= self.next_pos and self.next_pos > 0):
                self.update_features()
        if P < self.pos:
            while(self.cur_feat_pos and P < self.cur_feat_pos):
                self.update_features_backwards()
        self.pos = P
        self.seek_pos()
        pass

    def __del__(self):
        self.file.close()
        self.mt_file.close()

class View:
    def __init__(self, reader, stdscr):
        self.reader = reader
        self.screen = stdscr
        self.top_pos = self.reader.pos
        self.title_pos = 0

    def print_status(self):
        status = "{} ({:.3f}%)".format(self.top_pos, self.top_pos*100/self.reader.ch_size)
        if self.reader.current_info:
            status += " {} ({})".format(self.reader.current_info, strand_decode[self.reader.current_info_strand])

        self.screen.addstr(0, 0, status)

    def next_line(self):
        global scrx, scry, scrw, scrh, paused
        scrx = 0
        scry += 1
        if scry >= scrh:
            scry = 0
            return True
        else:
            return False

    def next_char(self):
        global scrx, scry, scrw, scrh, pos, size
        scrx += 1
        if scrx >= scrw-1:
            return self.next_line()
        else:
            return False

    def print_char(self, char, pair):
        global scrx, scry
        try:
            self.screen.addch(scry, scrx, char, curses.color_pair(pair))
            return self.next_char()
        except curses.error:
            return False
            #
            r = self.next_char()
            return r or self.print_char(char, pair)

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
        global prev_nucleotide, highlight
        pair = None
        nucleotide = self.reader.read()
        features = self.reader.current_features
        if feature_encode['gap'] in features:
            nucleotide = 4
            pair = PAIR_UNK
        elif feature_encode['CDS'] in features:
            if self.current_cds_phase is None:
                self.current_cds_phase = self.reader.get_cds_phase()
            if self.current_cds_phase & 4:
                pair = PAIR_CDS2 + nucleotide
            else:
                pair = PAIR_CDS + nucleotide
            self.current_cds_phase += 1
            if self.current_cds_phase & 3 == 3:
                self.current_cds_phase ^= 4
                self.current_cds_phase &= 4
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
        return self.print_char(nucleotide_decoding[nucleotide], pair)

    def print_title_line(self, title, line):
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
        height = len(chars['0'])
        n = height + line + 1
        if n >= 0 and n < height:
            length = 0
            for C in title:
                length += len(chars[C][0])
            pad = (scrw - length)//2
            if pad < 0:
                pad = 0
            for m in range(0, pad):
                self.print_char(" ", PAIR_UNK)
            for C in title:
                for c in chars[C][n]:
                    self.print_char(c, PAIR_UNK)
        else:
            for m in range(0, scrw):
                self.print_char(" ", PAIR_UNK)

    def fill(self):
        global scry
        self.current_cds_phase = None
        if self.title_pos < 0 and scry < -self.title_pos:
            self.print_title_line(self.reader.ch, self.title_pos + scry)
            return
        while not self.reader.eof:
            full = self.print_nucleotide()
            self.reader.advance()
            if full:
                break
        self.print_status()

    def scroll_down(self, n):
        global scrw, scrh, scrx, scry
        while n > 0 and not self.reader.eof:
            self.screen.scroll(1)
            if self.title_pos < 0:
                self.title_pos += 1
            if self.top_pos + (scrw-1)*(scrh+1) > 1:
                title = False
            else:
                title = True
            self.reader.jump_to(max(self.top_pos + (scrw-1)*scrh, 1))
            self.top_pos += scrw-1
            scry = scrh-1
            if self.top_pos + (scrw-1)*scrh < 1 and not title:
                scrx = 1 - self.top_pos
            else:
                scrx = 0
            self.fill()
            n -= 1
        if self.reader.current_info and self.reader.prev_info_pos and self.reader.prev_info_pos < self.top_pos:
            self.reader.current_info = ""

    def scroll_up(self, n):
        global scrw, scrh, scrx, scry
        while n > 0 and self.title_pos >= -10:
            self.screen.scroll(-1)
            tmp = scrh
            if self.top_pos <= 1:
                title = True
                if self.title_pos == 0:
                    self.reader.jump_to(max(self.top_pos, 1))
                    scrh = 1
                    scry = 1
                    scrx = 0
                    while scrx < 1 - self.top_pos:
                        self.print_char(' ', 0)
                    self.fill()
                    self.title_pos = -1
                else:
                    self.title_pos -= 1
            else:
                title = False
            self.top_pos -= scrw-1
            if not title:
                self.reader.jump_to(max(self.top_pos, 1))
            scrh = 2
            scry = 0
            if self.top_pos < 1 and not title:
                scrx = 1 - self.top_pos
            else:
                scrx = 0
            self.fill()
            scrh = tmp
            n -= 1
        if self.reader.current_info and self.reader.prev_info_pos and self.reader.prev_info_pos > self.top_pos + (scrw-1)*scrh:
            self.reader.current_info = ""

    def resize(self, W, H):
        global scrw, scrh, scrx, scry
        self.reader.jump_to(max(self.top_pos, 1))
        scrx = scry = 0
        scrw = W
        scrh = H
        self.fill()

    def __del__(self):
        pass

def main(stdscr):
    global paused, current_reader, scrw, scrh, ch_initial, pos_initial, pos_percent
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

    reader = Reader(ch_initial, pos_initial, pos_percent)
    view = View(reader, stdscr)
    view.fill()

    exit = False
    while not exit:
        if not paused:
            view.scroll_down(1)
            time.sleep(0.1)
        key = stdscr.getch()
        if not paused:
            while key in [curses.KEY_DOWN, curses.KEY_UP]:
                key = stdscr.getch()
        if key == curses.KEY_RESIZE:
            H, W = view.screen.getmaxyx()
            view.resize(W, H)
        elif key == ord('\n') or key == curses.KEY_ENTER or key == ord(' '):
            paused = not paused
        elif key == curses.KEY_DOWN:
            view.scroll_down(1)
        elif key == curses.KEY_UP:
            view.scroll_up(1)
        elif key == 27:
            exit = True

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
