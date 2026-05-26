"""
Real chromosome data for HaploTreeSim.

Chromosome lengths from hg38 (GRCh38) reference genome.
Source: UCSC Genome Browser
"""

# hg38 chromosome lengths in base pairs
# Source: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
HG38_CHROMOSOME_LENGTHS = {
    'chr1': 248956422,
    'chr2': 242193529,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
    'chr10': 133797422,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr13': 114364328,
    'chr14': 107043718,
    'chr15': 101991189,
    'chr16': 90338345,
    'chr17': 83257441,
    'chr18': 80373285,
    'chr19': 58617616,
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468,
    'chrX': 156040895,
    'chrY': 57227415,
}

AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
SEX_CHROMOSOMES = ['chrX', 'chrY']
WHOLE_GENOME_CHROMOSOMES = AUTOSOMES + SEX_CHROMOSOMES


def normalize_chromosomes(chromosomes) -> list:
    """Normalize a chromosome or chromosome list and validate hg38 names."""
    if isinstance(chromosomes, str):
        normalized = [chromosomes]
    else:
        normalized = list(chromosomes)

    if not normalized:
        raise ValueError("At least one chromosome must be specified")

    for chromosome in normalized:
        get_chromosome_length(chromosome)

    return normalized


def get_chromosome_length(chromosome: str) -> int:
    """Get the length of a chromosome in base pairs."""
    if chromosome not in HG38_CHROMOSOME_LENGTHS:
        raise ValueError(f"Unknown chromosome: {chromosome}")
    return HG38_CHROMOSOME_LENGTHS[chromosome]


def create_bins_for_chromosome(chromosome: str, bin_width: int) -> int:
    """Calculate number of bins needed for a chromosome."""
    length = get_chromosome_length(chromosome)
    num_bins = (length + bin_width - 1) // bin_width
    return num_bins


def create_bins_for_chromosomes(chromosomes, bin_width: int) -> int:
    """Calculate total bins needed for a chromosome list."""
    return sum(create_bins_for_chromosome(chromosome, bin_width) for chromosome in normalize_chromosomes(chromosomes))


def get_genome_length(chromosomes) -> int:
    """Return total reference length for a chromosome list."""
    return sum(get_chromosome_length(chromosome) for chromosome in normalize_chromosomes(chromosomes))


def get_arm_boundaries(chromosome: str) -> tuple:
    """Get approximate centromere position for chromosome arms."""
    CENTROMERES = {
        'chr1': (122026459, 125184587),
        'chr2': (92188145, 94090557),
        'chr3': (90772458, 93655574),
        'chr4': (49708101, 51743951),
        'chr5': (46485900, 50059807),
        'chr6': (58553888, 59829934),
        'chr7': (58169653, 61528020),
        'chr8': (43958052, 45338887),
        'chr9': (43389635, 45518558),
        'chr10': (39686682, 41593521),
        'chr11': (51078348, 54425074),
        'chr12': (34769407, 37185252),
        'chr13': (16000000, 18051248),
        'chr14': (16000000, 18173523),
        'chr15': (17083673, 19725254),
        'chr16': (36311158, 38265669),
        'chr17': (22813679, 26885980),
        'chr18': (15460899, 20861206),
        'chr19': (24498980, 27190874),
        'chr20': (26436232, 30038348),
        'chr21': (10864560, 12915808),
        'chr22': (12954789, 15054318),
        'chrX': (58605579, 62412542),
        'chrY': (10316944, 10544039),
    }
    
    if chromosome not in CENTROMERES:
        length = get_chromosome_length(chromosome)
        mid = length // 2
        return (mid - 1000000, mid + 1000000)
    
    return CENTROMERES[chromosome]


def describe_chromosome(chromosome: str, bin_width: int = 500000) -> dict:
    """Get detailed information about a chromosome."""
    length = get_chromosome_length(chromosome)
    num_bins = create_bins_for_chromosome(chromosome, bin_width)
    p_end, q_start = get_arm_boundaries(chromosome)
    
    p_arm_bins = p_end // bin_width
    q_arm_bins = (length - q_start) // bin_width
    
    return {
        'chromosome': chromosome,
        'length_bp': length,
        'length_mb': length / 1e6,
        'bin_width': bin_width,
        'num_bins': num_bins,
        'centromere_p_end': p_end,
        'centromere_q_start': q_start,
        'p_arm_bins': p_arm_bins,
        'q_arm_bins': q_arm_bins,
        'p_arm_length_mb': p_end / 1e6,
        'q_arm_length_mb': (length - q_start) / 1e6,
    }
