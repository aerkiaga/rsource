#!/usr/bin/python3

import os, sys, curses, time, configparser, re, bisect, shutil, copy, collections

chromosomes = [
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
    '11', '12', '13', '14', '15', '16', '17', '18', '19',
    '20', '21', '22', 'X', 'Y', 'mt'
]

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

# The pair numbers to which the different colors will be assigned
PAIR_UNK = 0
PAIR_HIGHLIGHT = 1
# All the following get 4 pairs each, one per nucleotide
PAIR_NONE = 8
PAIR_EXON_PSEUDO = 12
PAIR_UTR_GENE = 16
PAIR_CDS = 20
PAIR_CDS2 = 24
PAIR_INTRON = 28
PAIR_TRNA = 32
PAIR_RRNA = 36
PAIR_MIRNA = 40

#Foreground color
nucleotide_colors = {
    0 : 9,
    1 : 11,
    2 : 10,
    3 : 14,
    4 : 5
}

#Background color
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

#Other colors
other_colors = {
    PAIR_HIGHLIGHT : 11
}

script_path = os.path.realpath(__file__)
path = os.path.dirname(script_path)

scrw = 80
scrh = 25

highlight = {
    'cpg' : False,
    'tata' : False
}

#consensus sequence, max differences
highlighter = {
    'cpg' : ("CG", 0),
    'tata' : ("TATAWAWR", 1)
}

ch_initial = "1"
pos_initial = 1
pos_percent = False
paused = False

class Reader:
    ch_readers = {}

    @classmethod
    def get_ch_reader(cls, ch):
        if ch not in cls.ch_readers:
            return cls(ch, 0)
        return cls.ch_readers[ch]

    #apply feature (start/end of region) to self.current_features
    def apply_feature(self, feat):
        #end of region
        if feat & end_encode:
            if (feat & feature_mask) not in self.current_features:
                return #D
            self.current_features[feat & feature_mask] -= 1
            if self.current_features[feat & feature_mask] == 0:
                del self.current_features[feat & feature_mask]
        #start of region
        else:
            if (feat & feature_mask) not in self.current_features:
                self.current_features[feat & feature_mask] = 0
            self.current_features[feat & feature_mask] += 1

    #get extra info about the current feature in metadata file (gene name, CDS phase)
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

    #apply current feature, get its info, get next feature from metadata file
    def update_features(self):
        self.cur_feat_pos = self.next_pos
        self.apply_feature(self.next_feat)
        info = self.get_feature_info()
        if info:
            self.current_info_strand = ord(info[0])
            self.current_info = info[1:]
            self.prev_info_pos = self.next_pos
        dword = self.mt_file.read(4)
        if dword == b"":
            self.next_pos = None
            self.next_feat = None
        else:
            self.next_pos = int.from_bytes(dword, byteorder='little', signed=False)
            self.next_feat = int.from_bytes(self.mt_file.read(1), byteorder='little', signed=False)

    #seeks metadata file to previous feature
    #only seeks, no other side effect; returns type of feature it lands in
    def unget_feature(self):
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

    #updates self.current_features and seeks metadata file backwards
    #opposite of update_features()
    def update_features_backwards(self):
        self.next_pos = self.cur_feat_pos
        if self.next_feat is None:
            # we are at the end-of-file, after the new next's type
            self.mt_file.seek(self.mt_file.tell() - 1)
            lost_feat = int.from_bytes(self.mt_file.read(1), byteorder='little', signed=False) # new next (previous current)
        else:
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

    #updates self.current_features and seeks metadata to start
    def jump_to_mt_start(self):
        self.mt_file.seek(0)
        self.cur_feat_pos = None
        self.next_pos = int.from_bytes(self.mt_file.read(4), byteorder='little', signed=False)
        self.next_feat = int.from_bytes(self.mt_file.read(1), byteorder='little', signed=False)
        self.pos = 0
        self.current_features.clear()

    #updates self.current_features and seeks metadata to end
    def jump_to_mt_end(self):
        self.mt_file.seek(0, 2)
        self.next_pos = self.cur_feat_pos = None
        self.unget_feature()
        self.cur_feat_pos = int.from_bytes(self.mt_file.read(4), byteorder='little', signed=False)
        self.mt_file.seek(0, 2)
        self.next_feat = None
        self.current_features.clear()

    #gets CDS phase at self.pos
    def get_cds_phase(self):
        saved_fpos = self.mt_file.tell()
        #use cached CDS data if possible
        if saved_fpos == self.cds_phase_cache['saved_fpos']:
            start_phase = self.cds_phase_cache['start_phase']
            start_pos = self.cds_phase_cache['start_pos']
            relative_pos = self.pos - start_pos
            r = (start_phase + relative_pos%3) | (((relative_pos // 3) & 1) << 2)
        #find CDS and read its data, then restore file position
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

    #get byte from data file
    def get_byte(self):
        self.byte = self.file.read(1)
        if self.byte == b'':
            self.eof = True
        else:
            self.eof = False

    #seek to arbitrary chromosome position
    def seek_pos(self):
        self.file.seek(4 + (self.pos-1)//4)
        self.n = (self.pos-1) % 4
        self.get_byte()

    def __init__(self, ch, pos, pos_is_percent=False):
        self.ch = ch
        ch_path = os.path.join(path, self.ch + ".bin")
        self.file = open(ch_path, 'rb')
        ch_size = self.file.read(4)
        self.ch_size = int.from_bytes(ch_size, byteorder='little', signed=False)

        if pos_is_percent:
            pos = (pos * self.ch_size) // 100
        if pos_initial <= 0:
            pos = self.ch_size + pos + 1

        mt_path = os.path.join(path, self.ch + ".dat")
        self.mt_file = open(mt_path, 'rb')

        self.current_features = {}
        self.current_info = ""
        self.prev_info_pos = None

        self.jump_to_mt_start()
        self.jump_to(pos)

        self.cds_phase_cache = {
            'saved_fpos' : None,
            'start_phase' : None,
            'start_pos' : None
        }

        Reader.ch_readers[ch] = self

    #read nucleotide at current position
    def read(self):
        b = int.from_bytes(self.byte, byteorder='little')
        nucleotide = (b >> (2*(3-self.n))) & 3
        return nucleotide

    #advance to next nucleotide, without updating features
    def advance_nucleotide(self):
        self.n += 1
        if self.n & 3 == 0:
            self.get_byte()
            self.n &= 3
        self.last_nucleotides.append(self.read())
        self.pos += 1

    #advance to next nucleotide, possibly updating features
    def advance(self):
        global scrw, scrh
        self.advance_nucleotide()
        if self.pos == self.next_pos:
            while self.pos == self.next_pos:
                self.update_features()
        if self.pos > self.ch_size:
            self.eof = True

    #jump to arbitrary chromosome position, including any features
    def jump_to(self, P):
        if P == 1:
            self.jump_to_mt_start()
        elif P == self.ch_size:
            self.jump_to_mt_end()
        elif P < abs(P - self.pos):
            self.jump_to(1)
        elif (self.ch_size - P) < abs(P - self.pos):
            self.jump_to(self.ch_size)

        if P > self.pos:
            while(self.next_pos and P >= self.next_pos):
                self.update_features()
        if P < self.pos:
            while(self.cur_feat_pos and P < self.cur_feat_pos):
                self.update_features_backwards()
        self.last_nucleotides = collections.deque([None] * 20, 20)
        self.pos = P - 20 if P > 20 else 1
        self.seek_pos()
        while self.pos < P:
            self.advance_nucleotide()

    def __del__(self):
        self.file.close()
        self.mt_file.close()

class View:
    class Pos:
        def __init__(self):
            self.reader = None
            self.pos = None
            self.title_pos = None

        def __init__(self, reader, pos):
            self.reader = reader
            self.pos = pos
            self.title_pos = None

        #move reader to this position
        def sync_reader(self):
            if self.pos == self.reader.pos:
                pass
            elif self.pos == self.reader.pos+1:
                self.reader.advance()
            else:
                self.reader.jump_to(max(self.pos, 1))

        #get previous chromosome name
        def prev_ch_name(self):
            index = chromosomes.index(self.reader.ch)
            if index == 0:
                return None
            return chromosomes[index-1]

        #get next chromosome name
        def next_ch_name(self):
            index = chromosomes.index(self.reader.ch)
            if index == len(chromosomes)-1:
                return None
            return chromosomes[index+1]

        #check if viewer is in chromosome title
        def istitle(self):
            return self.title_pos is not None

        #check if viewer is in chromosome margin
        def ismargin(self):
            if self.pos < 1:
                return True
            self.sync_reader()
            return self.reader.eof

        #jump from end of chromosome to start of next
        def next_ch(self):
            self.pos = self.pos - self.reader.ch_size
            ch = self.next_ch_name()
            self.reader = Reader.get_ch_reader(ch)
            self.title_pos = -10
            self.sync_reader()

        #jump from start of chromosome to end of previous
        def prev_ch(self):
            ch = self.prev_ch_name()
            self.reader = Reader.get_ch_reader(ch)
            self.pos = self.reader.ch_size + self.pos
            if self.pos > self.reader.ch_size:
                self.pos -= scrw-1
            self.title_pos = None

        #advance view to next nucleotide
        #will not affect display unless we are at the margin of a chromosome
        def advance(self):
            self.pos += 1
            if not self.ismargin():
                if self.pos == 1:
                    self.reader.jump_to(1)
                else:
                    self.reader.advance()

        #move view to next line
        def next_line(self):
            if self.istitle():
                self.title_pos += 1
                if self.title_pos == 0:
                    self.title_pos = None
            else:
                if self.pos + scrw-1 > self.reader.ch_size and (self.next_ch_name() is not None):
                    self.next_ch()
                else:
                    self.pos += scrw-1

        #move view to previous line
        def prev_line(self):
            if self.pos <= 1:
                if not self.istitle():
                    self.title_pos = -1
                else:
                    if self.title_pos == -10:
                        self.prev_ch()
                    else:
                        self.title_pos -= 1
            else:
                self.pos -= scrw-1

        #advance a number of lines
        def advance_lines(self, n):
            for l in range(0, n):
                self.next_line()

        #check check if view should jump to next chromosome, and do it if so
        def check_ch_end(self):
            if self.pos > self.reader.ch_size and (self.next_ch_name() is not None):
                self.pos -= scrw-1
                self.next_ch()

        #check whether the view can be scrolled down
        def can_scroll_down(self):
            return (self.pos + (scrw-1)*scrh <= self.reader.ch_size) or (self.next_ch_name() is not None)

        #check whether the view can be scrolled up
        def can_scroll_up(self):
            return (not self.istitle()) or (self.title_pos > -10) or (self.prev_ch_name() is not None)

    def __init__(self, reader, stdscr):
        self.screen = stdscr
        self.top_pos = self.Pos(reader, reader.pos)

    #print status line on top
    def print_status(self):
        status = "{} ({:.3f}%)".format(self.top_pos.pos, self.top_pos.pos*100/self.top_pos.reader.ch_size)
        if self.top_pos.reader.current_info:
            status += " {} ({})".format(self.top_pos.reader.current_info, strand_decode[self.top_pos.reader.current_info_strand])

        self.screen.addstr(0, 0, status)

    #set the <number> characters before the last written one to a color pair
    def set_prev_pairs(self, number, pair):
        global scrw, scrh
        x = self.fillx
        y = self.filly
        for n in range(0, number):
            x -= 1
            if x < 0:
                x = scrw-2
                y -= 1
                if y < 0:
                    break
            self.screen.chgat(y, x, curses.color_pair(pair))

    #match nucleotide with consensus sequence symbol
    def match_consensus(self, nucleotide, consensus):
        nucleotide = nucleotide_decoding[nucleotide]
        return (
            nucleotide == consensus
            or consensus == 'N'
            or consensus == 'W' and nucleotide in ['A', 'T']
            or consensus == 'S' and nucleotide in ['C', 'G']
            or consensus == 'R' and nucleotide in ['A', 'G']
            or consensus == 'Y' and nucleotide in ['C', 'T']
            or consensus == 'M' and nucleotide in ['A', 'C']
            or consensus == 'K' and nucleotide in ['G', 'T']
            or consensus == 'B' and nucleotide in ['C', 'G', 'T']
            or consensus == 'D' and nucleotide in ['A', 'G', 'T']
            or consensus == 'H' and nucleotide in ['A', 'C', 'T']
            or consensus == 'V' and nucleotide in ['A', 'C', 'G']
        )

    #possibly apply highlight to preceding nucleotides, return pair for current one
    def apply_highlight(self, reader):
        global highlight, highlighter
        pair = None
        for name, enabled in highlight.items():
            if enabled:
                position = -1
                differences = 0
                while -position <= len(highlighter[name][0]):
                    if not self.match_consensus(reader.last_nucleotides[position], highlighter[name][0][position]):
                        differences += 1
                        if differences > highlighter[name][1]:
                            position = None
                            break
                    position -= 1
                if position is not None:
                    pair = PAIR_HIGHLIGHT
                    #also set the pair before this
                    self.set_prev_pairs(len(highlighter[name][0]) - 1, PAIR_HIGHLIGHT)
        return pair

    #get the appropriate nucleotide and pair for the current view position
    def get_nucleotide_and_pair(self, reader):
        pair = None
        nucleotide = reader.read()
        features = reader.current_features
        highlight_pair = self.apply_highlight(reader)
        if highlight_pair is not None:
            pair = highlight_pair
        else:
            if feature_encode['gap'] in features:
                nucleotide = 4
                pair = PAIR_UNK
            elif feature_encode['CDS'] in features:
                if self.current_cds_phase is None:
                    self.current_cds_phase = reader.get_cds_phase()
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
        return (nucleotide, pair)

    #print a line of the title of a chromosome
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
                self.fillx += 1
            for C in title:
                for c in chars[C][n]:
                    self.print_char(c, PAIR_UNK)
                    self.fillx += 1
            while self.fillx < scrw-1:
                self.print_char(" ", PAIR_UNK)
                self.fillx += 1
        else:
            for m in range(0, scrw-1):
                self.print_char(" ", PAIR_UNK)
                self.fillx += 1

    #write a character and pair to the current screen position
    def print_char(self, char, pair):
        try:
            self.screen.addch(self.filly, self.fillx, char, curses.color_pair(pair))
        except curses.error:
            return False

    #fill a portion of the screen with the appropriate characters and pairs
    def fill(self, x, y, h):
        global scrw
        self.fillx, self.filly = x, y
        self.fillmaxy = self.filly + h
        self.current_cds_phase = None

        #create a copy of the current view top position
        pos = copy.copy(self.top_pos)
        pos.advance_lines(y)

        while self.filly < self.fillmaxy:
            self.fillx = 0
            if pos.istitle():
                self.print_title_line(pos.reader.ch, pos.title_pos)
                pos.next_line()
            else:
                while self.fillx < scrw-1:
                    if pos.ismargin():
                        self.print_char(' ', 0)
                    else:
                        nucleotide, pair = self.get_nucleotide_and_pair(pos.reader)
                        self.print_char(nucleotide_decoding[nucleotide], pair)
                    pos.advance()
                    self.fillx += 1
                pos.check_ch_end()
            self.filly += 1
        self.print_status()

    #try to scroll view down a number of lines
    def scroll_down(self, n):
        global scrw, scrh
        while n > 0 and self.top_pos.can_scroll_down():
            self.screen.scroll(1)
            self.top_pos.next_line()
            self.fill(x=0, y=scrh-1, h=1)
            n -= 1
        if self.top_pos.reader.current_info and self.top_pos.reader.prev_info_pos and self.top_pos.reader.prev_info_pos < self.top_pos.pos:
            self.top_pos.reader.current_info = ""

    #try to scroll view up a number of lines
    def scroll_up(self, n):
        global scrw, scrh
        while n > 0 and self.top_pos.can_scroll_up():
            self.screen.scroll(-1)
            self.top_pos.prev_line()
            self.fill(x=0, y=0, h=2)
            n -= 1
        if self.top_pos.reader.current_info and self.top_pos.reader.prev_info_pos and self.top_pos.reader.prev_info_pos > self.top_pos.pos + (scrw-1)*scrh:
            self.top_pos.reader.current_info = ""

    #resize view to new size
    def resize(self, W, H):
        global scrw, scrh
        self.screen.clear()
        scrw = W
        scrh = H
        self.fill(x=0, y=0, h=scrh)

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
    view.fill(x=0, y=0, h=scrh)

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
