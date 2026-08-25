# ============================================================
# 🧬 GENELENS v2.0
# Neuroscience × Bioinformatics × Computational Biology
# ============================================================

import unittest
import math
from collections import Counter


# ============================================================
# 1. SEQUENCE PREPARATION
# ============================================================

def normalize_sequence(dna):
    # Remove spaces and convert the sequence to uppercase.
    return dna.strip().upper()


def sequence_length(dna):
    # Return the number of bases.
    return len(dna)


def count_bases(dna):
    # Count each standard DNA base.
    return {
        "A": dna.count("A"),
        "C": dna.count("C"),
        "G": dna.count("G"),
        "T": dna.count("T")
    }


def calculate_gc(dna):
    # Calculate GC percentage.
    if len(dna) == 0:
        return 0.0

    gc_count = dna.count("G") + dna.count("C")

    return (gc_count / len(dna)) * 100


def find_invalid_bases(dna):
    # Identify characters that are not standard DNA bases.
    valid_bases = {"A", "C", "G", "T"}

    invalid = []

    for base in dna:

        if base not in valid_bases and base not in invalid:
            invalid.append(base)

    return invalid


def analyze_sequence(dna):
    # Create a complete sequence profile.
    dna = normalize_sequence(dna)

    return {
        "sequence": dna,
        "length": sequence_length(dna),
        "counts": count_bases(dna),
        "gc_percent": calculate_gc(dna),
        "invalid_bases": find_invalid_bases(dna)
    }


# ============================================================
# 2. MOTIF ANALYSIS
# ============================================================

def find_motif(dna, motif):
    # Find every occurrence of a motif.
    dna = normalize_sequence(dna)
    motif = normalize_sequence(motif)

    positions = []

    if motif == "":
        return positions

    motif_length = len(motif)

    for i in range(len(dna) - motif_length + 1):

        if dna[i:i + motif_length] == motif:
            positions.append(i)

    return positions


# ============================================================
# 3. REVERSE COMPLEMENT
# ============================================================

def reverse_complement(dna):
    # Generate the reverse-complementary DNA strand.

    dna = normalize_sequence(dna)

    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    result = ""

    for base in dna[::-1]:

        if base not in complement:
            raise ValueError(
                f"Invalid DNA base: {base}"
            )

        result += complement[base]

    return result


# ============================================================
# 4. READING FRAME / STOP CODON ANALYSIS
# ============================================================

def find_stop_codons(dna, frame=0):
    # Find stop codons in a specific reading frame.

    dna = normalize_sequence(dna)

    stop_codons = {
        "TAA",
        "TAG",
        "TGA"
    }

    found = []

    for i in range(frame, len(dna) - 2, 3):

        codon = dna[i:i + 3]

        if codon in stop_codons:

            found.append({
                "codon": codon,
                "position": i
            })

    return found


# ============================================================
# 5. CODON TABLE
# ============================================================

CODON_TABLE = {

    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",

    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",

    "TAT": "Y",
    "TAC": "Y",

    "TAA": "*",
    "TAG": "*",
    "TGA": "*",

    "TGT": "C",
    "TGC": "C",
    "TGG": "W",

    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",

    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",

    "CAT": "H",
    "CAC": "H",

    "CAA": "Q",
    "CAG": "Q",

    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",

    "ATT": "I",
    "ATC": "I",
    "ATA": "I",

    "ATG": "M",

    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",

    "AAT": "N",
    "AAC": "N",

    "AAA": "K",
    "AAG": "K",

    "AGT": "S",
    "AGC": "S",

    "AGA": "R",
    "AGG": "R",

    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",

    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",

    "GAT": "D",
    "GAC": "D",

    "GAA": "E",
    "GAG": "E",

    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}


# ============================================================
# 6. DNA → PROTEIN TRANSLATION
# ============================================================

def translate_dna(dna, frame=0):
    # Translate DNA into amino-acid letters.

    dna = normalize_sequence(dna)

    protein = ""

    for i in range(frame, len(dna) - 2, 3):

        codon = dna[i:i + 3]

        if codon in CODON_TABLE:
            protein += CODON_TABLE[codon]

        else:
            protein += "X"

    return protein


# ============================================================
# 7. OPEN READING FRAME ANALYSIS
# ============================================================

def find_orfs(dna, minimum_length=30):
    # Find ORFs beginning with ATG and ending at a stop codon.

    dna = normalize_sequence(dna)

    orfs = []

    for frame in range(3):

        i = frame

        while i <= len(dna) - 3:

            codon = dna[i:i + 3]

            if codon == "ATG":

                start = i
                j = i + 3

                while j <= len(dna) - 3:

                    stop_codon = dna[j:j + 3]

                    if stop_codon in {
                        "TAA",
                        "TAG",
                        "TGA"
                    }:

                        length = j + 3 - start

                        if length >= minimum_length:

                            sequence = dna[
                                start:j + 3
                            ]

                            orfs.append({
                                "frame": frame,
                                "start": start,
                                "end": j + 3,
                                "length": length,
                                "stop": stop_codon,
                                "sequence": sequence,
                                "protein": translate_dna(
                                    sequence,
                                    0
                                )
                            })

                        i = j + 3
                        break

                    j += 3

                else:
                    i += 3

            else:
                i += 3

    return orfs


# ============================================================
# 8. K-MER ANALYSIS
# ============================================================

def kmer_counts(dna, k):
    # Count all k-length subsequences.

    dna = normalize_sequence(dna)

    if k <= 0:
        return {}

    if k > len(dna):
        return {}

    kmers = []

    for i in range(len(dna) - k + 1):

        kmers.append(
            dna[i:i + k]
        )

    return dict(
        Counter(kmers)
    )


# ============================================================
# 9. SLIDING-WINDOW GC ANALYSIS
# ============================================================

def sliding_gc(dna, window=20):
    # Calculate GC content across overlapping windows.

    dna = normalize_sequence(dna)

    if window <= 0:
        return []

    if window > len(dna):
        return []

    results = []

    for start in range(
        len(dna) - window + 1
    ):

        section = dna[
            start:start + window
        ]

        gc = (
            section.count("G")
            +
            section.count("C")
        )

        gc_percent = (
            gc / window
        ) * 100

        results.append({
            "start": start,
            "end": start + window,
            "gc": gc_percent
        })

    return results


# ============================================================
# 10. MUTATION CONSEQUENCE ANALYSIS
# ============================================================

def classify_mutation(
    original,
    mutated,
    position,
    frame=0
):
    # Classify a single-base substitution.

    original = normalize_sequence(original)
    mutated = normalize_sequence(mutated)

    if len(original) != len(mutated):

        return {
            "type": "length_change"
        }

    if position < 0 or position >= len(original):

        return {
            "type": "invalid_position"
        }

    if original[position] == mutated[position]:

        return {
            "type": "no_change"
        }

    original_protein = translate_dna(
        original,
        frame
    )

    mutated_protein = translate_dna(
        mutated,
        frame
    )

    codon_position = (
        position - frame
    ) // 3

    if codon_position < 0:

        return {
            "type": "outside_frame"
        }

    if codon_position >= len(
        original_protein
    ):

        return {
            "type": "outside_frame"
        }

    original_aa = original_protein[
        codon_position
    ]

    mutated_aa = mutated_protein[
        codon_position
    ]

    if mutated_aa == "*":

        mutation_type = "nonsense"

    elif original_aa == mutated_aa:

        mutation_type = "synonymous"

    else:

        mutation_type = "missense"

    return {
        "type": mutation_type,
        "position": position,
        "original_base": original[position],
        "mutated_base": mutated[position],
        "original_amino_acid": original_aa,
        "mutated_amino_acid": mutated_aa
    }


# ============================================================
# 11. FASTA READER
# ============================================================

def read_fasta(filename):
    # Read sequences from a multi-FASTA file.

    sequences = {}

    current_id = None
    current_sequence = ""

    try:

        with open(
            filename,
            "r"
        ) as file:

            for line in file:

                line = line.strip()

                if line == "":
                    continue

                if line.startswith(">"):

                    if current_id is not None:

                        sequences[
                            current_id
                        ] = current_sequence

                    header = line[1:]

                    current_id = (
                        header.split()[0]
                    )

                    current_sequence = ""

                else:

                    current_sequence += (
                        line.upper()
                    )

            if current_id is not None:

                sequences[
                    current_id
                ] = current_sequence

    except FileNotFoundError:

        print(
            "⚠️ FASTA file was not found."
        )

        return {}

    except Exception as error:

        print(
            "⚠️ Could not read FASTA file."
        )

        print(
            "Error:",
            error
        )

        return {}

    return sequences


# ============================================================
# 12. NEUROSCIENCE GENE KNOWLEDGE BASE
# ============================================================

NEUROSCIENCE_GENES = {

    "SYN1": {
        "category": "Synaptic",
        "system": "Neural communication",
        "role": "Synaptic vesicle-associated protein"
    },

    "BDNF": {
        "category": "Neuroplasticity",
        "system": "Learning and memory",
        "role": "Neurotrophin involved in neuronal signaling"
    },

    "SCN1A": {
        "category": "Neuronal",
        "system": "Electrical signaling",
        "role": "Voltage-gated sodium channel"
    },

    "GRIN1": {
        "category": "Synaptic",
        "system": "Neural communication",
        "role": "NMDA receptor subunit"
    },

    "MAPT": {
        "category": "Neurodegeneration",
        "system": "Cytoskeletal regulation",
        "role": "Tau protein"
    },

    "APP": {
        "category": "Neurodegeneration",
        "system": "Neuronal biology",
        "role": "Amyloid precursor protein"
    },

    "RAD51": {
        "category": "DNA repair",
        "system": "Genome stability",
        "role": "Homologous recombination"
    },

    "RPA1": {
        "category": "DNA repair",
        "system": "Genome stability",
        "role": "Replication protein A"
    },

    "SRSF1": {
        "category": "RNA biology",
        "system": "RNA processing",
        "role": "RNA splicing factor"
    }
}


# ============================================================
# 13. BASIC SEQUENCE ANALYSIS
# ============================================================

def full_sequence_analysis():

    print()
    print("🧬 FULL SEQUENCE ANALYSIS")
    print("-" * 50)

    dna = input(
        "Enter your DNA sequence: "
    )

    results = analyze_sequence(dna)

    if results["length"] == 0:

        print(
            "⚠️ No sequence entered."
        )

        return

    print()
    print("=" * 60)
    print("🧬 GENELENS SEQUENCE PROFILE")
    print("=" * 60)

    print()
    print("SEQUENCE")
    print("-" * 30)
    print(results["sequence"])

    print()
    print("📏 LENGTH")
    print("-" * 30)
    print(results["length"], "bp")

    print()
    print("🧬 BASE COMPOSITION")
    print("-" * 30)

    counts = results["counts"]

    print("A:", counts["A"])
    print("C:", counts["C"])
    print("G:", counts["G"])
    print("T:", counts["T"])

    print()
    print("📊 GC CONTENT")
    print("-" * 30)

    print(
        round(
            results["gc_percent"],
            2
        ),
        "%"
    )

    print()
    print("⚠️ VALIDATION")
    print("-" * 30)

    if len(results["invalid_bases"]) == 0:

        print(
            "✓ No invalid bases detected."
        )

    else:

        print(
            "Invalid bases:",
            results["invalid_bases"]
        )

        return

    print()
    print("🔄 REVERSE COMPLEMENT")
    print("-" * 30)

    print(
        reverse_complement(
            results["sequence"]
        )
    )

    print()
    print("🛑 READING FRAME ANALYSIS")
    print("-" * 30)

    for frame in range(3):

        stops = find_stop_codons(
            results["sequence"],
            frame
        )

        print()
        print(
            "Reading frame",
            frame
        )

        if len(stops) == 0:

            print(
                "No stop codons found."
            )

        else:

            for stop in stops:

                print(
                    "Stop:",
                    stop["codon"],
                    "| Position:",
                    stop["position"]
                )

    print()
    print("=" * 60)


# ============================================================
# 14. MOTIF ANALYSIS
# ============================================================

def motif_analysis():

    print()
    print("🔎 MOTIF ANALYSIS")
    print("-" * 50)

    dna = input(
        "Enter your DNA sequence: "
    )

    dna = normalize_sequence(dna)

    if dna == "":

        print(
            "⚠️ No sequence entered."
        )

        return

    invalid = find_invalid_bases(dna)

    if len(invalid) > 0:

        print(
            "⚠️ Invalid bases:",
            invalid
        )

        return

    motif = input(
        "Enter motif to search for: "
    )

    motif = normalize_sequence(motif)

    if motif == "":

        print(
            "⚠️ No motif entered."
        )

        return

    positions = find_motif(
        dna,
        motif
    )

    print()
    print(
        "Motif:",
        motif
    )

    print(
        "Occurrences:",
        len(positions)
    )

    if positions:

        print(
            "Positions:",
            positions
        )

    else:

        print(
            "Motif was not found."
        )


# ============================================================
# 15. REVERSE COMPLEMENT ANALYSIS
# ============================================================

def reverse_complement_analysis():

    print()
    print("🔄 REVERSE COMPLEMENT")
    print("-" * 50)

    dna = normalize_sequence(
        input(
            "Enter your DNA sequence: "
        )
    )

    if dna == "":

        print(
            "⚠️ No sequence entered."
        )

        return

    invalid = find_invalid_bases(dna)

    if len(invalid) > 0:

        print(
            "⚠️ Invalid bases:",
            invalid
        )

        return

    print()
    print("Original:")
    print(dna)

    print()
    print("Reverse complement:")
    print(
        reverse_complement(dna)
    )


# ============================================================
# 16. STOP CODON ANALYSIS
# ============================================================

def stop_codon_analysis():

    print()
    print("🛑 READING FRAME ANALYSIS")
    print("-" * 50)

    dna = normalize_sequence(
        input(
            "Enter your DNA sequence: "
        )
    )

    if dna == "":
        print("⚠️ No sequence entered.")
        return

    invalid = find_invalid_bases(dna)

    if len(invalid) > 0:

        print(
            "⚠️ Invalid bases:",
            invalid
        )

        return

    for frame in range(3):

        stops = find_stop_codons(
            dna,
            frame
        )

        print()
        print(
            "Reading frame:",
            frame
        )

        if not stops:

            print(
                "No stop codons found."
            )

        else:

            for stop in stops:

                print(
                    stop["codon"],
                    "at position",
                    stop["position"]
                )


# ============================================================
# 17. SEQUENCE INFORMATION
# ============================================================

def sequence_information():

    print()
    print("📋 SEQUENCE INFORMATION")
    print("-" * 50)

    dna = input(
        "Enter your DNA sequence: "
    )

    results = analyze_sequence(dna)

    if results["length"] == 0:

        print(
            "⚠️ No sequence entered."
        )

        return

    print()
    print("Sequence:", results["sequence"])
    print("Length:", results["length"], "bp")
    print(
        "GC:",
        round(
            results["gc_percent"],
            2
        ),
        "%"
    )

    print(
        "Invalid bases:",
        results["invalid_bases"]
    )


# ============================================================
# 18. NEUROSCIENCE GENE ANALYSIS
# ============================================================

def neuroscience_gene_analysis():

    print()
    print("🧠 NEUROSCIENCE GENE ANALYSIS")
    print("-" * 55)

    gene = input(
        "Enter gene symbol: "
    ).strip().upper()

    if gene == "":

        print(
            "⚠️ No gene symbol entered."
        )

        return

    print()
    print("Choose biological category:")
    print("1. 🧠 Neuronal")
    print("2. 🔗 Synaptic")
    print("3. 🌱 Neurodevelopment")
    print("4. 🧬 Other")

    category_choice = input(
        "Choose category (1-4): "
    ).strip()

    categories = {
        "1": "Neuronal",
        "2": "Synaptic",
        "3": "Neurodevelopment",
        "4": "Other"
    }

    if category_choice not in categories:

        print(
            "⚠️ Invalid category."
        )

        return

    dna = normalize_sequence(
        input(
            "Enter reference DNA sequence: "
        )
    )

    if dna == "":

        print(
            "⚠️ No sequence entered."
        )

        return

    results = analyze_sequence(dna)

    if results["invalid_bases"]:

        print(
            "⚠️ Invalid bases:",
            results["invalid_bases"]
        )

        return

    print()
    print("=" * 65)
    print("🧠 GENELENS NEUROSCIENCE GENE PROFILE")
    print("=" * 65)

    print()
    print("Gene:", gene)
    print(
        "Category:",
        categories[category_choice]
    )

    print()
    print("Length:", results["length"], "bp")

    print(
        "GC content:",
        round(
            results["gc_percent"],
            2
        ),
        "%"
    )

    print()
    print("Base composition:")

    for base, count in results[
        "counts"
    ].items():

        print(
            base + ":",
            count
        )

    print()
    print("Protein translation, frame 0:")

    print(
        translate_dna(
            dna,
            0
        )
    )

    print()
    print("ORFs:")

    orfs = find_orfs(
        dna,
        30
    )

    print(
        "Detected:",
        len(orfs)
    )

    print()
    print("=" * 65)


# ============================================================
# 19. NEUROSCIENCE EXPRESSION ANALYSIS
# ============================================================

def neuroscience_expression_analysis():

    print()
    print("🧠📊 NEUROSCIENCE EXPRESSION ANALYSIS")
    print("-" * 55)

    gene = input(
        "Enter gene symbol: "
    ).strip().upper()

    if gene == "":

        print(
            "⚠️ No gene symbol entered."
        )

        return

    regions = [
        "Cortex",
        "Cerebellum",
        "Hippocampus",
        "Brainstem"
    ]

    expression = {}

    for region in regions:

        try:

            expression[region] = float(
                input(
                    f"{region} expression: "
                )
            )

        except ValueError:

            print(
                "⚠️ Values must be numbers."
            )

            return

    highest_region = max(
        expression,
        key=expression.get
    )

    lowest_region = min(
        expression,
        key=expression.get
    )

    mean_expression = (
        sum(expression.values())
        /
        len(expression)
    )

    print()
    print("=" * 60)
    print("🧠📊 EXPRESSION PROFILE")
    print("=" * 60)

    print()
    print("Gene:", gene)

    for region, value in expression.items():

        print(
            region + ":",
            value
        )

    print()
    print(
        "Highest:",
        highest_region,
        expression[highest_region]
    )

    print(
        "Lowest:",
        lowest_region,
        expression[lowest_region]
    )

    print(
        "Mean:",
        round(
            mean_expression,
            3
        )
    )

    print(
        "Range:",
        round(
            expression[highest_region]
            -
            expression[lowest_region],
            3
        )
    )

    print()
    print("=" * 60)


# ============================================================
# 20. COMPARE NEUROSCIENCE GENES
# ============================================================

def compare_neuroscience_genes():

    print()
    print("🔬 COMPARE NEUROSCIENCE GENES")
    print("-" * 55)

    try:

        number = int(
            input(
                "How many genes? "
            )
        )

    except ValueError:

        print(
            "⚠️ Enter a whole number."
        )

        return

    if number < 2:

        print(
            "⚠️ At least 2 genes are required."
        )

        return

    genes = []

    for i in range(number):

        print()
        print(
            "GENE",
            i + 1
        )

        gene = input(
            "Gene symbol: "
        ).strip().upper()

        dna = normalize_sequence(
            input(
                "DNA sequence: "
            )
        )

        if dna == "":

            print(
                "⚠️ Sequence cannot be empty."
            )

            return

        results = analyze_sequence(dna)

        if results["invalid_bases"]:

            print(
                "⚠️ Invalid bases:",
                results["invalid_bases"]
            )

            return

        genes.append({
            "gene": gene,
            "length": results["length"],
            "gc": results["gc_percent"]
        })

    print()
    print("=" * 65)
    print("🔬 GENE COMPARISON")
    print("=" * 65)

    print()

    for item in genes:

        print(
            f"{item['gene']:<15}"
            f"{item['length']:<12}"
            f"{item['gc']:.2f}%"
        )

    longest = max(
        genes,
        key=lambda x: x["length"]
    )

    highest_gc = max(
        genes,
        key=lambda x: x["gc"]
    )

    print()
    print(
        "Longest:",
        longest["gene"],
        longest["length"],
        "bp"
    )

    print(
        "Highest GC:",
        highest_gc["gene"],
        round(
            highest_gc["gc"],
            2
        ),
        "%"
    )

    print()
    print("=" * 65)


# ============================================================
# 21. FASTA DATASET ANALYSIS
# ============================================================

def fasta_dataset_analysis():

    print()
    print("🧬 FASTA DATASET ANALYSIS")
    print("-" * 55)

    filename = input(
        "FASTA file path: "
    ).strip()

    sequences = read_fasta(
        filename
    )

    if not sequences:

        print(
            "⚠️ No sequences found."
        )

        return

    results = []

    for sequence_id, dna in sequences.items():

        analysis = analyze_sequence(dna)

        results.append({
            "id": sequence_id,
            "length": analysis["length"],
            "gc": analysis["gc_percent"],
            "invalid": analysis["invalid_bases"]
        })

    print()
    print("=" * 75)
    print("🧬 GENELENS DATASET PROFILE")
    print("=" * 75)

    print()

    print(
        f"{'ID':<20}"
        f"{'Length':<12}"
        f"{'GC %':<12}"
        f"{'Status':<15}"
    )

    print("-" * 75)

    for result in results:

        status = (
            "Valid"
            if not result["invalid"]
            else "Invalid"
        )

        print(
            f"{result['id']:<20}"
            f"{result['length']:<12}"
            f"{result['gc']:<12.2f}"
            f"{status:<15}"
        )

    valid = [
        x for x in results
        if not x["invalid"]
    ]

    if not valid:

        print(
            "⚠️ No valid sequences."
        )

        return

    average_length = (
        sum(
            x["length"]
            for x in valid
        )
        /
        len(valid)
    )

    average_gc = (
        sum(
            x["gc"]
            for x in valid
        )
        /
        len(valid)
    )

    longest = max(
        valid,
        key=lambda x: x["length"]
    )

    highest_gc = max(
        valid,
        key=lambda x: x["gc"]
    )

    print()
    print("📊 DATASET SUMMARY")
    print("-" * 45)

    print(
        "Total sequences:",
        len(results)
    )

    print(
        "Valid sequences:",
        len(valid)
    )

    print(
        "Average length:",
        round(
            average_length,
            2
        ),
        "bp"
    )

    print(
        "Average GC:",
        round(
            average_gc,
            2
        ),
        "%"
    )

    print(
        "Longest:",
        longest["id"],
        longest["length"],
        "bp"
    )

    print(
        "Highest GC:",
        highest_gc["id"],
        round(
            highest_gc["gc"],
            2
        ),
        "%"
    )


# ============================================================
# 22. DATASET VISUALIZATION
# ============================================================

def visualize_neuroscience_dataset():

    print()
    print("📈 DATASET VISUALIZATION")
    print("-" * 55)

    filename = input(
        "FASTA file path: "
    ).strip()

    sequences = read_fasta(
        filename
    )

    if not sequences:

        print(
            "⚠️ No sequences found."
        )

        return

    results = []

    for sequence_id, dna in sequences.items():

        analysis = analyze_sequence(dna)

        if analysis["invalid_bases"]:
            continue

        results.append({
            "id": sequence_id,
            "length": analysis["length"],
            "gc": analysis["gc_percent"],
            "counts": analysis["counts"]
        })

    if not results:

        print(
            "⚠️ No valid sequences."
        )

        return

    while True:

        print()
        print("📈 VISUALIZATION MENU")
        print("-" * 45)

        print(
            "1. 📊 GC-content comparison"
        )

        print(
            "2. 📏 Sequence-length comparison"
        )

        print(
            "3. 🧬 Base composition"
        )

        print(
            "4. 📋 Dataset summary"
        )

        print(
            "5. ↩️ Back"
        )

        choice = input(
            "Choose (1-5): "
        ).strip()

        if choice == "1":

            print()
            print("📊 GC CONTENT")
            print("-" * 55)

            for result in results:

                bar = "█" * int(
                    result["gc"] / 2
                )

                print(
                    f"{result['id']:<20}"
                    f"{result['gc']:>7.2f}% "
                    f"{bar}"
                )

        elif choice == "2":

            print()
            print("📏 SEQUENCE LENGTH")
            print("-" * 55)

            maximum = max(
                x["length"]
                for x in results
            )

            for result in results:

                bar_length = int(
                    (
                        result["length"]
                        /
                        maximum
                    )
                    * 40
                )

                print(
                    f"{result['id']:<20}"
                    f"{result['length']:>8} bp "
                    f"{'█' * bar_length}"
                )

        elif choice == "3":

            print()
            print("🧬 BASE COMPOSITION")
            print("-" * 55)

            for result in results:

                print()
                print(
                    result["id"]
                )

                for base, count in result[
                    "counts"
                ].items():

                    print(
                        f"  {base}: {count}"
                    )

        elif choice == "4":

            average_gc = (
                sum(
                    x["gc"]
                    for x in results
                )
                /
                len(results)
            )

            average_length = (
                sum(
                    x["length"]
                    for x in results
                )
                /
                len(results)
            )

            print()
            print("📋 DATASET SUMMARY")
            print("-" * 45)

            print(
                "Sequences:",
                len(results)
            )

            print(
                "Average GC:",
                round(
                    average_gc,
                    2
                ),
                "%"
            )

            print(
                "Average length:",
                round(
                    average_length,
                    2
                ),
                "bp"
            )

        elif choice == "5":

            break

        else:

            print(
                "⚠️ Invalid option."
            )


# ============================================================
# 23. ORF + PROTEIN ANALYSIS
# ============================================================

def orf_analysis():

    print()
    print("🧬 ORF + PROTEIN ANALYSIS")
    print("-" * 55)

    dna = normalize_sequence(
        input(
            "Enter DNA sequence: "
        )
    )

    if dna == "":

        print(
            "⚠️ No sequence entered."
        )

        return

    invalid = find_invalid_bases(dna)

    if invalid:

        print(
            "⚠️ Invalid bases:",
            invalid
        )

        return

    minimum_text = input(
        "Minimum ORF length in bp (default 30): "
    ).strip()

    try:

        minimum_length = (
            int(minimum_text)
            if minimum_text
            else 30
        )

    except ValueError:

        print(
            "⚠️ Invalid minimum length."
        )

        return

    orfs = find_orfs(
        dna,
        minimum_length
    )

    print()
    print("=" * 70)
    print("🧬 GENELENS ORF REPORT")
    print("=" * 70)

    print(
        "Sequence length:",
        len(dna),
        "bp"
    )

    print(
        "ORFs detected:",
        len(orfs)
    )

    if not orfs:

        print(
            "No ORFs meeting the criteria were detected."
        )

        return

    for number, orf in enumerate(
        orfs,
        1
    ):

        print()
        print(
            "ORF",
            number
        )

        print("-" * 45)

        print(
            "Frame:",
            orf["frame"]
        )

        print(
            "Start:",
            orf["start"]
        )

        print(
            "End:",
            orf["end"]
        )

        print(
            "Length:",
            orf["length"],
            "bp"
        )

        print(
            "Stop:",
            orf["stop"]
        )

        print(
            "DNA:",
            orf["sequence"]
        )

        print(
            "Protein:",
            orf["protein"]
        )


# ============================================================
# 24. K-MER ANALYSIS
# ============================================================

def kmer_analysis():

    print()
    print("🔢 K-MER ANALYSIS")
    print("-" * 55)

    dna = normalize_sequence(
        input(
            "Enter DNA sequence: "
        )
    )

    if dna == "":
        print("⚠️ No sequence entered.")
        return

    if find_invalid_bases(dna):

        print(
            "⚠️ Invalid DNA sequence."
        )

        return

    try:

        k = int(
            input(
                "Enter k (example: 3): "
            )
        )

    except ValueError:

        print(
            "⚠️ k must be an integer."
        )

        return

    counts = kmer_counts(
        dna,
        k
    )

    if not counts:

        print(
            "⚠️ No k-mers available."
        )

        return

    print()
    print("=" * 60)
    print("🔢 K-MER FREQUENCY PROFILE")
    print("=" * 60)

    sorted_kmers = sorted(
        counts.items(),
        key=lambda item: item[1],
        reverse=True
    )

    for kmer, count in sorted_kmers:

        print(
            f"{kmer:<10}{count}"
        )

    print()
    print(
        "Unique k-mers:",
        len(counts)
    )

    print(
        "Total k-mers:",
        sum(counts.values())
    )


# ============================================================
# 25. SLIDING-WINDOW GC ANALYSIS
# ============================================================

def sliding_gc_analysis():

    print()
    print("📊 SLIDING-WINDOW GC ANALYSIS")
    print("-" * 55)

    dna = normalize_sequence(
        input(
            "Enter DNA sequence: "
        )
    )

    if dna == "":
        print("⚠️ No sequence entered.")
        return

    try:

        window = int(
            input(
                "Window size (default 20): "
            ) or 20
        )

    except ValueError:

        print(
            "⚠️ Window must be an integer."
        )

        return

    results = sliding_gc(
        dna,
        window
    )

    if not results:

        print(
            "⚠️ Window is larger than sequence."
        )

        return

    highest = max(
        results,
        key=lambda x: x["gc"]
    )

    lowest = min(
        results,
        key=lambda x: x["gc"]
    )

    average = (
        sum(
            x["gc"]
            for x in results
        )
        /
        len(results)
    )

    print()
    print("=" * 60)
    print("📊 GC LANDSCAPE")
    print("=" * 60)

    print(
        "Windows analyzed:",
        len(results)
    )

    print(
        "Average GC:",
        round(
            average,
            2
        ),
        "%"
    )

    print(
        "Highest GC:",
        round(
            highest["gc"],
            2
        ),
        "%",
        "at position",
        highest["start"]
    )

    print(
        "Lowest GC:",
        round(
            lowest["gc"],
            2
        ),
        "%",
        "at position",
        lowest["start"]
    )

    print()

    for result in results[:50]:

        bar = "█" * int(
            result["gc"] / 2
        )

        print(
            f"{result['start']:>5} "
            f"{result['gc']:>6.2f}% "
            f"{bar}"
        )


# ============================================================
# 26. MUTATION CONSEQUENCE SIMULATOR
# ============================================================

def mutation_analysis():

    print()
    print("🧪 MUTATION CONSEQUENCE SIMULATOR")
    print("-" * 60)

    original = normalize_sequence(
        input(
            "Original DNA sequence: "
        )
    )

    if original == "":
        print("⚠️ No sequence entered.")
        return

    if find_invalid_bases(original):

        print(
            "⚠️ Invalid DNA sequence."
        )

        return

    try:

        position = int(
            input(
                "Position to mutate (0-based): "
            )
        )

    except ValueError:

        print(
            "⚠️ Position must be an integer."
        )

        return

    if position < 0 or position >= len(original):

        print(
            "⚠️ Position outside sequence."
        )

        return

    new_base = input(
        "New base (A/C/G/T): "
    ).strip().upper()

    if new_base not in {
        "A",
        "C",
        "G",
        "T"
    }:

        print(
            "⚠️ Invalid base."
        )

        return

    mutated = (
        original[:position]
        +
        new_base
        +
        original[position + 1:]
    )

    result = classify_mutation(
        original,
        mutated,
        position
    )

    print()
    print("=" * 65)
    print("🧪 MUTATION REPORT")
    print("=" * 65)

    print(
        "Original base:",
        result.get(
            "original_base",
            "N/A"
        )
    )

    print(
        "New base:",
        result.get(
            "mutated_base",
            "N/A"
        )
    )

    print(
        "Position:",
        result.get(
            "position",
            "N/A"
        )
    )

    print()
    print(
        "Mutation type:",
        result["type"]
    )

    if "original_amino_acid" in result:

        print(
            "Original amino acid:",
            result["original_amino_acid"]
        )

        print(
            "New amino acid:",
            result["mutated_amino_acid"]
        )

    print()
    print("Original:")
    print(original)

    print()
    print("Mutated:")
    print(mutated)

    print()
    print("=" * 65)


# ============================================================
# 27. BAYESIAN NEURAL DECODER
# ============================================================

def gaussian_probability(
    x,
    mean,
    standard_deviation
):
    # Gaussian likelihood function.

    if standard_deviation <= 0:

        return 0.0

    coefficient = (
        1
        /
        (
            standard_deviation
            *
            math.sqrt(
                2 * math.pi
            )
        )
    )

    exponent = (
        -(
            (x - mean) ** 2
        )
        /
        (
            2
            *
            standard_deviation ** 2
        )
    )

    return (
        coefficient
        *
        math.exp(exponent)
    )


def bayesian_decoder():

    print()
    print("🧠 BAYESIAN NEURAL DECODER")
    print("-" * 60)

    print()
    print(
        "Estimate which stimulus is most likely"
    )

    print(
        "given an observed neural response."
    )

    print()

    try:

        response = float(
            input(
                "Observed neural response: "
            )
        )

        mean_a = float(
            input(
                "Mean response for stimulus A: "
            )
        )

        sd_a = float(
            input(
                "Standard deviation for A: "
            )
        )

        mean_b = float(
            input(
                "Mean response for stimulus B: "
            )
        )

        sd_b = float(
            input(
                "Standard deviation for B: "
            )
        )

        prior_a = float(
            input(
                "Prior probability of A: "
            )
        )

        prior_b = float(
            input(
                "Prior probability of B: "
            )
        )

    except ValueError:

        print(
            "⚠️ All values must be numerical."
        )

        return

    if prior_a < 0 or prior_b < 0:

        print(
            "⚠️ Prior probabilities cannot be negative."
        )

        return

    likelihood_a = gaussian_probability(
        response,
        mean_a,
        sd_a
    )

    likelihood_b = gaussian_probability(
        response,
        mean_b,
        sd_b
    )

    weighted_a = (
        likelihood_a
        *
        prior_a
    )

    weighted_b = (
        likelihood_b
        *
        prior_b
    )

    evidence = (
        weighted_a
        +
        weighted_b
    )

    if evidence == 0:

        print(
            "⚠️ Posterior cannot be calculated."
        )

        return

    posterior_a = (
        weighted_a
        /
        evidence
    )

    posterior_b = (
        weighted_b
        /
        evidence
    )

    print()
    print("=" * 70)
    print("🧠 BAYESIAN DECODING RESULT")
    print("=" * 70)

    print()
    print(
        "Observed response:",
        response
    )

    print(
        "P(response | A):",
        round(
            likelihood_a,
            6
        )
    )

    print(
        "P(response | B):",
        round(
            likelihood_b,
            6
        )
    )

    print()
    print(
        "P(A | response):",
        round(
            posterior_a,
            4
        )
    )

    print(
        "P(B | response):",
        round(
            posterior_b,
            4
        )
    )

    print()

    if posterior_a > posterior_b:

        print(
            "🧠 Prediction: STIMULUS A"
        )

    elif posterior_b > posterior_a:

        print(
            "🧠 Prediction: STIMULUS B"
        )

    else:

        print(
            "🧠 Prediction: TIE"
        )

    print()
    print("=" * 70)


# ============================================================
# 28. GENE KNOWLEDGE EXPLORER
# ============================================================

def gene_knowledge_explorer():

    print()
    print("🧠 GENE KNOWLEDGE EXPLORER")
    print("-" * 55)

    gene = input(
        "Enter gene symbol: "
    ).strip().upper()

    if gene not in NEUROSCIENCE_GENES:

        print()
        print(
            "⚠️ Gene is not currently in the GeneLens knowledge base."
        )

        print()
        print(
            "Available genes:"
        )

        print(
            ", ".join(
                NEUROSCIENCE_GENES.keys()
            )
        )

        return

    data = NEUROSCIENCE_GENES[gene]

    print()
    print("=" * 65)
    print("🧬 GENELENS BIOLOGICAL PROFILE")
    print("=" * 65)

    print()
    print(
        "Gene:",
        gene
    )

    print(
        "Category:",
        data["category"]
    )

    print(
        "System:",
        data["system"]
    )

    print(
        "Role:",
        data["role"]
    )

    print()
    print("=" * 65)


# ============================================================
# 29. RESEARCH PROJECT DASHBOARD
# ============================================================

def research_dashboard():

    print()
    print("=" * 75)
    print("🧬 GENELENS RESEARCH PROJECT")
    print("   Neuroscience × Bioinformatics × Computation")
    print("=" * 75)

    print()
    print("PROJECT QUESTION")
    print("-" * 55)

    print(
        "How can computational methods help us"
    )

    print(
        "understand biological systems and neural function?"
    )

    print()
    print("LEARNING → APPLICATION")
    print("-" * 55)

    print(
        "🧠 Neuroscience"
    )

    print(
        "   Neural communication, sensory systems,"
    )

    print(
        "   plasticity, neurodevelopment and decision-making."
    )

    print()

    print(
        "🧬 Bioinformatics"
    )

    print(
        "   DNA sequences, FASTA, motifs, ORFs,"
    )

    print(
        "   proteins, mutations and k-mers."
    )

    print()

    print(
        "💻 Python"
    )

    print(
        "   Functions, data structures, file processing"
    )

    print(
        "   and computational analysis."
    )

    print()

    print(
        "📊 Probability"
    )

    print(
        "   Bayesian inference and neural decoding."
    )

    print()

    print(
        "🔬 Research"
    )

    print(
        "   Question → Data → Analysis → Interpretation."
    )

    print()
    print("CURRENT GENELENS CAPABILITIES")
    print("-" * 55)

    capabilities = [
        "DNA sequence analysis",
        "GC-content analysis",
        "Motif detection",
        "Reverse complements",
        "Reading-frame analysis",
        "DNA → protein translation",
        "Open reading frame detection",
        "k-mer analysis",
        "Sliding-window GC analysis",
        "Mutation consequence simulation",
        "FASTA dataset analysis",
        "Neuroscience gene profiles",
        "Expression comparison",
        "Bayesian neural decoding",
        "Dataset visualization",
        "Automated unit testing"
    ]

    for capability in capabilities:

        print(
            "✓",
            capability
        )

    print()
    print("=" * 75)


# ============================================================
# 30. UNIT TESTS
# ============================================================

class TestGeneLens(unittest.TestCase):

    def test_sequence_length(self):

        self.assertEqual(
            sequence_length("ATGC"),
            4
        )

    def test_base_counts(self):

        self.assertEqual(
            count_bases("AATCGG"),
            {
                "A": 2,
                "C": 1,
                "G": 2,
                "T": 1
            }
        )

    def test_gc_content(self):

        self.assertAlmostEqual(
            calculate_gc("ATGC"),
            50.0
        )

    def test_empty_gc(self):

        self.assertEqual(
            calculate_gc(""),
            0.0
        )

    def test_invalid_bases(self):

        self.assertEqual(
            find_invalid_bases("ATGX"),
            ["X"]
        )

    def test_reverse_complement(self):

        self.assertEqual(
            reverse_complement("ATGC"),
            "GCAT"
        )

    def test_motif(self):

        self.assertEqual(
            find_motif(
                "ATGATG",
                "ATG"
            ),
            [0, 3]
        )

    def test_translation(self):

        self.assertEqual(
            translate_dna("ATGGCC"),
            "MA"
        )

    def test_stop_translation(self):

        self.assertEqual(
            translate_dna("ATGTAA"),
            "M*"
        )

    def test_kmers(self):

        self.assertEqual(
            kmer_counts(
                "ATAT",
                2
            ),
            {
                "AT": 2,
                "TA": 1
            }
        )

    def test_sliding_gc(self):

        results = sliding_gc(
            "GGCCAT",
            4
        )

        self.assertAlmostEqual(
            results[0]["gc"],
            100.0
        )

    def test_orf(self):

        orfs = find_orfs(
            "ATGAAATAA",
            9
        )

        self.assertEqual(
            len(orfs),
            1
        )

    def test_mutation(self):

        original = "ATGGCTTAA"
        mutated = "ATGACTTAA"

        result = classify_mutation(
            original,
            mutated,
            4
        )

        self.assertEqual(
            result["type"],
            "missense"
        )

    def test_gaussian_probability(self):

        result = gaussian_probability(
            0,
            0,
            1
        )

        self.assertGreater(
            result,
            0
        )


# ============================================================
# 31. MAIN MENU
# ============================================================

def main():

    while True:

        print()

        print(
            "╔════════════════════════════════════════════════════╗"
        )

        print(
            "║                 🧬 GENELENS                        ║"
        )

        print(
            "║       NEUROSCIENCE × BIOINFORMATICS                ║"
        )

        print(
            "║     Place the discover the beautiful world of      ║"
        )
        print(
            "║         neuroscience and bioinformatics            ║"
        )

        print(
            "╚════════════════════════════════════════════════════╝"
        )

        print()

        print("MAIN MENU")
        print("-" * 60)

        print(
            "1.  📊 Full sequence analysis"
        )

        print(
            "2.  🔎 Search for a motif"
        )

        print(
            "3.  🔄 Reverse complement"
        )

        print(
            "4.  🛑 Reading-frame / stop codons"
        )

        print(
            "5.  🧠 Neuroscience gene analysis"
        )

        print(
            "6.  🧠 Neuroscience expression analysis"
        )

        print(
            "7.  🔬 Compare neuroscience genes"
        )

        print(
            "8.  🧬 FASTA dataset analysis"
        )

        print(
            "9.  📈 Visualize neuroscience dataset"
        )

        print(
            "10. 📋 Sequence information"
        )

        print()
        print("─── GENELENS RESEARCH TOOLS ───")

        print(
            "11. 🧬 ORF + protein analysis"
        )

        print(
            "12. 🔢 k-mer analysis"
        )

        print(
            "13. 📊 Sliding-window GC analysis"
        )

        print(
            "14. 🧪 Mutation consequence simulator"
        )

        print(
            "15. 🧠 Bayesian neural decoder"
        )

        print(
            "16. 🧠 Gene knowledge explorer"
        )

        print(
            "17. 🌷 Research project dashboard"
        )

        print(
            "18. 🧪 Run tests"
        )

        print(
            "19. 🚪 Exit"
        )

        print()

        choice = input(
            "Choose an option (1-19): "
        ).strip()

        # ----------------------------------------------------
        # BASIC ANALYSIS
        # ----------------------------------------------------

        if choice == "1":

            full_sequence_analysis()

        elif choice == "2":

            motif_analysis()

        elif choice == "3":

            reverse_complement_analysis()

        elif choice == "4":

            stop_codon_analysis()

        elif choice == "5":

            neuroscience_gene_analysis()

        elif choice == "6":

            neuroscience_expression_analysis()

        elif choice == "7":

            compare_neuroscience_genes()

        elif choice == "8":

            fasta_dataset_analysis()

        elif choice == "9":

            visualize_neuroscience_dataset()

        elif choice == "10":

            sequence_information()

        # ----------------------------------------------------
        # RESEARCH TOOLS
        # ----------------------------------------------------

        elif choice == "11":

            orf_analysis()

        elif choice == "12":

            kmer_analysis()

        elif choice == "13":

            sliding_gc_analysis()

        elif choice == "14":

            mutation_analysis()

        elif choice == "15":

            bayesian_decoder()

        elif choice == "16":

            gene_knowledge_explorer()

        elif choice == "17":

            research_dashboard()

        # ----------------------------------------------------
        # TESTS
        # ----------------------------------------------------

        elif choice == "18":

            print()
            print(
                "🧪 RUNNING GENELENS TESTS"
            )

            print(
                "-" * 45
            )

            unittest.main(
                argv=[
                    "first-arg-is-ignored"
                ],
                exit=False
            )

        # ----------------------------------------------------
        # EXIT
        # ----------------------------------------------------

        elif choice == "19":

            print()
            print(
                "🧬 Thank you for using GeneLens."
            )

            print(
                "Keep asking biological questions. 💙"
            )

            break

        else:

            print()
            print(
                "⚠️ Invalid option."
            )

            print(
                "Please choose a number from 1 to 19."
            )


# ============================================================
# 32. RUN GENELENS
# ============================================================

if __name__ == "__main__":

    main()
