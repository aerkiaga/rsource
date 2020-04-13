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
landmarks anywhere. A negative value will count from the end of the chromosome.

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
