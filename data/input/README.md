# Data provenance

All data originates from GERMS Lab Box. Paths below are structured for syncing with `rclone`. The examples assume you have a remote named `Box:` configured in `rclone config`.

Optional convenience variables:

``` bash
BOX_BASE="Box:GERMS-DATAMAN/DOE-CABBI"
LOCAL_BASE="/home/baponte/gdrive_local/post_doc/DOE-CABBI/bolivar/germs_miscanthus/data/input"
```

## Summary

| Dataset | Origin (Box) | Destination (local) |
|------------------------|------------------------|------------------------|
| LAMPS — RNA amplicon sequencing | Box:GERMS-DATAMAN/DOE-CABBI/LAMPS/RNA\\ amplicon\\ sequencing/ | /home/baponte/gdrive_local/post_doc/DOE-CABBI/bolivar/germs_miscanthus/data/input/LAMPS/RNA_amplicon_sequencing/ |
| LAMPS — DNA amplicon sequencing | Box:GERMS-DATAMAN/DOE-CABBI/LAMPS/DNA\\ amplicon\\ sequencing/ | /home/baponte/gdrive_local/post_doc/DOE-CABBI/bolivar/germs_miscanthus/data/input/LAMPS/DNA_amplicon_sequencing/ |
| Energy Farm Collab — 16S | Box:GERMS-DATAMAN/DOE-CABBI/Energy Farm Collab/R/files_for_phyloseq_16S/ | /home/baponte/gdrive_local/post_doc/DOE-CABBI/bolivar/germs_miscanthus/data/input/energy_farm_collab/files_for_phyloseq_16S/ |
| Energy Farm Collab — AMF | Box:GERMS-DATAMAN/DOE-CABBI/Energy Farm Collab/R/files_for_phyloseq_AMF/ | /home/baponte/gdrive_local/post_doc/DOE-CABBI/bolivar/germs_miscanthus/data/input/energy_farm_collab/files_for_phyloseq_AMF/ |

Notes: - Corrected destination typos for clarity and consistency: - LAMPS DNA destination: added missing slash and fixed spelling to `.../LAMPS/DNA_amplicon_sequencing/` (was `LAMPSDNA_amplicon_sequnencing`). - Energy Farm AMF destination: fixed to `.../files_for_phyloseq_AMF/` (was `files_for_phylose_AMF`).

## Details and `rclone` commands

Tip: Use quotes to avoid escaping spaces in Box paths.

### LAMPS data

-   RNA amplicon sequencing
    -   Origin: `Box:GERMS-DATAMAN/DOE-CABBI/LAMPS/RNA\ amplicon\ sequencing/`

    -   Destination: `/home/baponte/gdrive_local/post_doc/DOE-CABBI/bolivar/germs_miscanthus/data/input/LAMPS/RNA_amplicon_sequencing/`

    -   Sync:

        ``` bash
        rclone sync "$BOX_BASE/LAMPS/RNA amplicon sequencing/" \
          "$LOCAL_BASE/LAMPS/RNA_amplicon_sequencing/" \
          --progress --checksum
        ```
-   DNA amplicon sequencing
    -   Origin: `Box:GERMS-DATAMAN/DOE-CABBI/LAMPS/DNA amplicon sequencing/`

    -   Destination: `/home/baponte/gdrive_local/post_doc/DOE-CABBI/bolivar/germs_miscanthus/data/input/LAMPS/DNA_amplicon_sequencing/`

    -   Sync:

        ``` bash
        rclone sync "$BOX_BASE/LAMPS/DNA amplicon sequencing/" \
          "$LOCAL_BASE/LAMPS/DNA_amplicon_sequencing/" \
          --progress --checksum
        ```

### Energy Farm Collab

-   16S
    -   Origin: `Box:GERMS-DATAMAN/DOE-CABBI/Energy Farm Collab/R/files_for_phyloseq_16S/`

    -   Destination: `/home/baponte/gdrive_local/post_doc/DOE-CABBI/bolivar/germs_miscanthus/data/input/energy_farm_collab/files_for_phyloseq_16S/`

    -   Sync:

        ``` bash
        rclone sync "$BOX_BASE/Energy Farm Collab/R/files_for_phyloseq_16S/" \
          "$LOCAL_BASE/energy_farm_collab/files_for_phyloseq_16S/" \
          --progress --checksum
        ```
-   AMF
    -   Origin: `Box:GERMS-DATAMAN/DOE-CABBI/Energy Farm Collab/R/files_for_phyloseq_AMF/`

    -   Destination: `/home/baponte/gdrive_local/post_doc/DOE-CABBI/bolivar/germs_miscanthus/data/input/energy_farm_collab/files_for_phyloseq_AMF/`

    -   Sync:

        ``` bash
        rclone sync "$BOX_BASE/Energy Farm Collab/R/files_for_phyloseq_AMF/" \
          "$LOCAL_BASE/energy_farm_collab/files_for_phyloseq_AMF/" \
          --progress --checksum
        ```

## Conventions

-   Trailing slash on directories is intentional for `rclone` to copy contents into the destination folder.
-   Local destination names use underscores for readability; Box paths retain original names.
-   Use `--dry-run` with `rclone` first to preview changes.
