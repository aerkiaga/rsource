# rsource

This is a console tool that downloads and displays the entire human genome,
highlighting various features throughout the 3+ billion base-pairs that comprise
the extremely obfuscated source code of our own species. It can also be customized
through a simple configuration file, and used as a didactic tool, a way to gain
insight into the genome, or just an overly complex console toy.

## Installing
The Python scripts can be run in place. They have the following dependencies:
 * Python 3 (3.5 or above is recommended)
 * PyPy3 (optional, but speeds up setup if detected)
 * wget
 * gunzip
 * rm (optional)

Most Linux systems already have these installed (save for PyPy3).
See [Configuration](#configuration) for setup customization. In order to set up the tool, one must run:

    python3 ./rsource.py

This will download **1 GB** of data, then process it into about **800 MB**, and
delete the original files. All these files will be created in the same directory
where the Python scripts are. Then, the script will exit. This step will be
performed only once, as long as the new files aren't deleted or renamed. This is
discouraged, as the setup process is long and requires a ridiculously large download
from the National Center for Biotechnology Information.

800 MB is the approximate size of the genome: 3 billion base-pairs (bp), where
a base-pair represents 2 bits of information. Our cells store all this data in
several DNA fragments called chromosomes, named 1 to 22, plus X and (optionally) Y;
we also have DNA in our mitochondria (mtDNA). The program stores the packed
sequences in files ending in `.bin`, and metadata (for highlighting) in files
ending in `.dat`. Since the sequence used (Genome Reference Consortium Human Build 38)
contains gaps, these are annotated in `.gap` files, then included in `.dat` files.
The `.gap` files can be safely deleted afterwards, although they don't take up
much space (about 7 kB).

## Running
After the setup step has been completed, the same script can be run as:

    python3 ./rsource.py

This will open the viewer, and start enumerating from the tip of the short arm
of chromosome 1, along the (+)-sense strand (our DNA has two "complementary" strands).
After a 10000 bp gap (unknown DNA), known regions are reached. As the status bar
in the top-left corner shows, chromosome 1 is huge. It would take days to reach
the end at typical speed, and more than a month to read all 25 files.

Fortunately, it's possible to have a look at any chromosome like this:

    python3 ./rsource.py 18

Possible values are 1 to 22, X, Y and mt. Furthermore, to start at a particular
location in the chromosome, one can use:

    python3 ./rsource.py 18.10000
    python3 ./rsource.py 18.50%
    python3 ./rsource.py 18.-700000
    python3 ./rsource.py 18.-1%

This advances a number of base-pairs into the chromosome (starts at 1), or a
percentage of its length, allowing to easily see particularly interesting
landmarks anywhere (see [Travel Guide](#travel-guide)). A negative value will
count from the end of the chromosome instead.

## Configuration
All configuration is done via a single `config.ini` file. This contains a few
sections with different options. If any option (or the whole file) is missing or
commented out with '#' at the start of a line, a default value will be used instead.

The "**Setup**" region contains settings that affect the setup step:
 * *delete sequence*: delete the downloaded sequence file (~950 MB) after setup.
 * *delete annotations*: delete the downloaded annotations file (~50 MB) after setup.
 * *delete gaps*: delete all generated `.gap` files (~7 kB) after setup.

The "**Nucleobase Colors**" section can be used to set foreground colors for the
different nucleobases: A (adenine), C (cytosine), G (guanine) and T (thymine).
Colors can be in HTML HEX or RGB format.

The "**Region Colors**" section sets background colors for highlighting different
regions. Colors can be in HTML HEX or RGB format.

## Travel Guide
Since our genome is so large, it's important to know where to search for interesting
items. Here is a list of regions to have a look at, laid out as a tutorial.

 * [Telomeres](#telomeres)
 * [Centromeres](#centromeres)
 * [Genes](#genes)
 * [Gene desert](#gene-desert)

### Telomeres
Run:

    python3 ./rsource.py 5.10000

Soon after the start, there is a long string of repeated bases. The repeated
sequence is `CCCTAA`. This sequence probably also extends into the "unknown" area.
Now run:

    python3 ./rsource.py 1.-11000

There are also repeated sequences before the large gap at the end. This time the
repeated fragment is `TTAGGG`.

There is a reason for this. Our DNA has two strands, and each of them has two ends,
termed 5' and 3'. The sense we are reading it is 5' to 3' (also called (+) sense).
The two strands are joined in opposite senses, with their bases paired A to T and C to G.
When DNA is copied, a new complementary strand is made from each of the two, resulting
in two double-stranded DNA helices. The machinery used for this is about the same
in all organisms, including bacteria (which have circular chromosomes with no ends)
and us. However, due to a quirk in it, a few bases in the two 3' ends of our
linear chromosomes are not copied, so the 5' end of the new complementary strand
is shorter.

This poses two problems. Firstly, after many copying operations, important regions
could be lost. Secondly, since only one strand is shortened at each end, non-paired
DNA results; this could pair randomly with similar sequences in other chromosomes
and cause trouble.

The way these problems are solved is by having a special protein (telomerase)
add these fixed sequences at the ends of all chromosomes. This "useless" DNA can
be lost without problem. Additionally, other proteins recognize these particular
repeated sequences, bind them, and "tie" them in a "hairpin" shape that leaves
no unpaired DNA exposed.

In all vertebrates, including humans, the repeated sequence is `TTAGGG`. That is,
reading from the center towards the ends, along the (+) strand in each end.
This is the sequence we find at the end of chromosome 1, since we are reading in
that sense. At the beginning of chromosome 5, however, we are reading *towards*
the center, so the sequence is reversed. Additionally, since our (+) strand is
the (-) strand counting from the center, we also get the complementary sequence.
The result is `CCCTAA`.

### Centromeres
Resize your console the nearest you can to 170 width. If you can't, resize it
as close as possible to 85 width (actually, any multiple or divisor of 170 will do,
if it is close enough and reasonably large). Make it as large as possible in height.
Now run:

    python3 ./rsource.py 1.50%

The resulting pattern is truly beautiful, and, more importantly, it's real DNA
inside our cells! The region of chromosome 1 where we have landed is called the
centromere. It just happens to be at the center of the chromosome, which is not
necessarily true for all chromosomes.

When a cell wants to divide, it packs all its DNA (except mtDNA) tightly, making
each chromosome into a rod-shaped bundle (these are the funky bars that appeared
during setup, with bright and dark bands that are seen when staining and observing
through a microscope). Since the cell has usually copied all chromosomes before
doing this, the two copies of each packed chromosome (called chromatids) are joined
together at the centromere, forming the characteristic 'X'-shaped chromosome.
Then, the cell pulls each chromatid from its centromere, separating them and
moving each to one pole. This way, chromosomes are sorted, one copy of each for
each of the two daughter cells.

What we see at this location is several thousand copies of a ~170 bp sequence
called an α-satellite. This is recognized by cell machinery to identify where
the centromere is. Note that copies of the α-satellite are slightly different,
both in length and content. Further along the centromere, lines on the console
start to become wobblier: a few α-satellites of different lengths are grouped, and
these groups themselves are also periodically repeated! Everything here has appeared
through mutation and copying, all throughout evolution, introducing randomness and
sheer complexity even in something as simple as a repeated sequence.

### Genes
Run:

    python3 ./rsource.py 16.172000

You should be right before ("upstream" of) the *HBA2* gene, for "Hemoglobin subunit
alpha 2". You may see the *HBA1* gene ("Hemoglobin subunit alpha 1") a few kbp
downstream. These two genes contain instructions to make these two subunits of
the larger protein hemoglobin. One α1, one α2 and two β subunits (the HBB gene
is in chromosome 11) form a complete hemoglobin molecule. Around these two genes
there are other hemoglobin genes: downstream we find *HBQ1* (θ1), upstream there
are *HBAP1* (a pseudogene, inactive remnant of a once real gene), *HBM* (μ), *HBZP1*
(another pseudogene) and *HBZ* (ζ, is part of embryo hemoglobin). The fact that
similar genes are right after one another is not a coincidence; they were produced
through duplications and mutations, forming what is known as a gene cluster.

Roughly speaking, genes contain instructions for making proteins, which are the
actual "workers" of the cell. Just as DNA is a string of nucleotides, proteins
are strings of amino acids, or assemblies of multiple such strings called subunits.
The process for making proteins from genes is really complex and differs between
organisms, but is fundamentally laid out like this: the DNA is *transcribed* into
a similar molecule termed RNA, this RNA (primary transcript) is edited to make
a messenger RNA (mRNA), and this mRNA is then read in chunks of 3 nucleotides
called "codons", and for each possible codon a corresponding amino acid (there are
20 of them) is appended to the growing amino acid chain.

Look at the gene. It seems to be made of three large chunks, separated by two
regions of a different type. When the gene is transcribed, the processing that
makes that primary transcript into a mRNA involves (among other things) the removal
of these intermediate regions (*introns*) and joing together of the three relevant
fragments (*exons*). This is called "splicing" and means that a large portion of
genetic material doesn't even get translated into protein. It serves a function,
though. In some cases (not this), the same primary transcript can be spliced in
different ways depending on the situation, leaving out a different set of regions
to get a different mRNA, and thus a different protein. It is generally observed
that more complex organisms tend to make more extensive use of this alternative
splicing.

After some other modifications that are not relevant to understanding the structure
of a gene, an mRNA is obtained. One important difference with RNA in general vs DNA
is that it has U (uracil) in place of T, with the same function. When the mRNA is
to be read to make a protein, the translation machinery looks for a "start codon".
This codon is AUG, so look for an ATG near the start. It's not *right* at the start,
leaving an untraslated region of mRNA (the 5' UTR). The start codon also codes for
an amino acid, the first of every translated protein. Then, 141 further codons are
encountered, making a 142 amino acid chain. Translation ends when one of the three
possible "stop codons" is found; in this gene, this last codon is UAA (TAA here).
The stop codon doen't code for any amino acid, it only signals the end of the
translated region and releases the new protein. After it, we can see another
untranslated sequence (the 3' UTR).

How does the transcription machinery know where to begin reading? If you look
closely a little upstream of the gene, you should find the sequence `CAAT`. This
is part of the longer sequence `CAGCCAATGA`, which signals a transcription site.
The sequence, called a "CAAT box", changes a lot between genes, but the `CAAT`
part is highly conserved. Other sequences decide how often and in what situations
transcription should start; all this elements form the *promoter* of a gene, and
are recognized by different proteins that control transcription.

There are also sequences to signal intron/exon boundaries. All the exons here that
are followed by an intron end with `AG`. After this, the intron starts with the
sequence `GTGAG`. Then, near the end of the intron, there is a sequence like
`CTCAC` or `CTGAC`, followed by a short run of mostly Cs and Ts, and finishing
with `CAG`. These sequences are not exactly the same in all genes (although very
similar), and there are alternate methods for delimiting introns and exons. The
basic idea, however, is usually the same: the ends of an intron are recognized by
proteins, brought together, and cut from the surrounding exons, that are joined.
Not that splicing is done after transcription, so Ts here would actually be Us.

The 5' UTR contains special sequences too. In fact, the start codon is part of a
larger sequence, `GAACCCACC ATG G` in this case, named Kozak sequence after its
discoverer, Marilyn S. Kozak. This ensures that the correct start location has
been chosen, even if other start codons are present in the mRNA.

Finally, there can also be special elements in the 3' UTR. Near the end we see the
sequence `AATAAA`, which basically signals that the transcript should be edited
to protect it from being destroyed right away.

### Gene desert
Run:

    python3 ./rsource.py 8.77500000

The presented sequence is located between two genes (*PEX2* and *PKIA*) that are
very far apart from each other. In fact, over a million bp away. From here, there
are hundreds of thousands of bp to the next feature.

What is its function? According to current knowledge, none. These regions comprise
about 72% of our genome, and not only most of them haven't been found to be useful
in any way, but mutations in these rarely appear to cause any effect whatsoever.
Along with an additional 26% of non-coding DNA *within* genes (introns), this
adds up to a ridiculous 98% of total potentially useless material.

Not all of this is "junk", however. Some of this non-coding DNA is known to perform
various functions. A number of heritable diseases result from mutations outside
genes, in these areas previously though to be irrelevant. Furthermore, a fraction
of non-coding material is conserved between species, which means that evolution
doesn't allow mutations in it "for free", suggesting some functional role.

Still, the total estimated "useful" DNA in our cells is just about 5-15% of the total.
The present region, for example, or at least most of it, is likely unnecessary.

Not all organisms have as much non-coding DNA as us (98%). Bacteria, for instance,
typically only have about 1/5 of it in their genome, in contrast to most other
organisms. Some, however, have apparently managed to remove it altogether; the
plant *Utricularia gibba* only has about 3%!

Have a look at the sequence. What probably stands out more in it is the presence
of relatively long strings of a single base, or a couple of them repeated over
and over. These are called microsatellites, and give us a glimpse of what our DNA
is actually made of.

We talk about an issue with DNA copying machinery in [Telomeres](#telomeres).
Here is another; when the copying proteins come across a series of short repeated
sequences, the strand being built can slip and pair with the wrong repeat. If this
happens, the resulting copy will have one repeat more or less than the original;
generally, the tendency is to increase the repeat number. This process is overall
highly unlikely, yet extremely common compared to other kinds of mutations. The
mutation rate further increases with repeat number, logarithmically.

This exemplifies an important truth about our DNA: most of it isn't there because
we need it, but because it is (or once was) somehow able to copy itself around
over and over. Microsatellites are a 3% of our genome, by themselves more than
coding DNA; the rest is no different, composed of sequences with various strategies
to "reproduce" inside the genome. DNA regions, thus, evolve by itself, for their
own purposes, as long as they don't cause much harm to the organism they exist in.
This is the Selfish Gene Theory of Natural Selection.
