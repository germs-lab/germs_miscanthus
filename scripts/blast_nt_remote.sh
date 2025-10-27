#!/usr/bin/env bash
set -Eeuo pipefail


# Usage: 
# ./blast_nt_remote.sh ASV.fasta blast_remote.tsv


# Remote BLAST against nt with an Entrez query excluding uninformative titles.
REGION="$1"
FASTA="$2"
OUT="$3"
#THREADS="${3:-2}"
FIELDS="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids ssciname stitle"


if [[ -z "$REGION" || -z "$FASTA" || -z "$OUT" ]]; then
  echo "Usage: $0 <REGION> <asv.fasta> <blast_output.tsv>" >&2
  echo "REGION: 16S, ITS, or AMF" >&2
  exit 1
fi

# Define region-specific parameters
if [[ "$REGION" == "16S" ]]; then
  ENTREZ_QUERY='"Bacteria"[Organism] NOT (uncultured[All Fields] OR environmental[All Fields] OR metagenom*[All Fields] OR clone[All Fields] OR unidentified[All Fields])'
  EVALUE="1e-10"
  PERC_IDENTITY="97"
  QCOV_HSP_PERC="85"
  MAX_TARGET_SEQS="50"
  
elif [[ "$REGION" == "ITS" ]]; then
  ENTREZ_QUERY='"fungi"[Organism] NOT (uncultured[All Fields] OR environmental[All Fields] OR metagenom*[All Fields] OR clone[All Fields] OR unidentified[All Fields])'
  EVALUE="1e-8"
  PERC_IDENTITY="95"
  QCOV_HSP_PERC="75"
  MAX_TARGET_SEQS="75"
  
elif [[ "$REGION" == "AMF" ]]; then
  ENTREZ_QUERY='"Glomeromycota"[Organism] NOT (uncultured[All Fields] OR environmental[All Fields] OR metagenom*[All Fields] OR clone[All Fields] OR unidentified[All Fields])'
  EVALUE="1e-6"
  PERC_IDENTITY="90"
  QCOV_HSP_PERC="70"
  MAX_TARGET_SEQS="100"
  
else
  echo "Error: REGION must be 16S, ITS, or AMF" >&2
  exit 1
fi

# Run BLAST with region-specific parameters
blastn \
  -task megablast \
  -db nt \
  -remote \
  -query "$FASTA" \
  -evalue "$EVALUE" \
  -perc_identity "$PERC_IDENTITY" \
  -qcov_hsp_perc "$QCOV_HSP_PERC" \
  -max_target_seqs "$MAX_TARGET_SEQS" \
  -max_hsps 1 \
  -dust no \
  -entrez_query "$ENTREZ_QUERY" \
  -outfmt "6 $FIELDS" \
  -out "$OUT"

tmp=$(mktemp)
printf '%s\n' "$(echo "$FIELDS" | tr ' ' '\t')" > "$tmp"
cat "$OUT" >> "$tmp"
mv "$tmp" "$OUT"



# BLAST 2.13.0+
# USAGE
#   blastn [-h] [-help] [-import_search_strategy filename]
#     [-export_search_strategy filename] [-task task_name] [-db database_name]
#     [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
#     [-negative_gilist filename] [-negative_seqidlist filename]
#     [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
#     [-negative_taxidlist filename] [-entrez_query entrez_query]
#     [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
#     [-subject subject_input_file] [-subject_loc range] [-query input_file]
#     [-out output_file] [-evalue evalue] [-word_size int_value]
#     [-gapopen open_penalty] [-gapextend extend_penalty]
#     [-perc_identity float_value] [-qcov_hsp_perc float_value]
#     [-max_hsps int_value] [-xdrop_ungap float_value] [-xdrop_gap float_value]
#     [-xdrop_gap_final float_value] [-searchsp int_value]
#     [-sum_stats bool_value] [-penalty penalty] [-reward reward] [-no_greedy]
#     [-min_raw_gapped_score int_value] [-template_type type]
#     [-template_length int_value] [-dust DUST_options]
#     [-filtering_db filtering_database]
#     [-window_masker_taxid window_masker_taxid]
#     [-window_masker_db window_masker_db] [-soft_masking soft_masking]
#     [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
#     [-best_hit_score_edge float_value] [-subject_besthit]
#     [-window_size int_value] [-off_diagonal_range int_value]
#     [-use_index boolean] [-index_name string] [-lcase_masking]
#     [-query_loc range] [-strand strand] [-parse_deflines] [-outfmt format]
#     [-show_gis] [-num_descriptions int_value] [-num_alignments int_value]
#     [-line_length line_length] [-html] [-sorthits sort_hits]
#     [-sorthsps sort_hsps] [-max_target_seqs num_sequences]
#     [-num_threads int_value] [-mt_mode int_value] [-remote] [-version]

