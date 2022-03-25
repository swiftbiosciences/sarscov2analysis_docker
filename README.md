# Swift Sarscov2 analysis script

## Pulling image from Docker hub
``` sh
docker pull swiftbiosci/sarscov2analysis:latest
```

## Download Masterfile from Box

``` sh
# Swift Normalase Amplicon SARS-CoV-2 (SNAP) With additinoal genome coverage Masterfile
https://danaherlifesciences.ent.box.com/s/lxzgbnfoc9zsycjhd647ijbxh5rad0m3
```

## Clone repository swiftbiosciences/sarscov2analysis_docker:

``` sh
git clone https://github.com/swiftbiosciences/sarscov2analysis_docker.git
```

## Wrapper script usage:
1. Copy run_swift_sarscov2_docker.sh to location of your working PATH.
2. Place all fastq files and panel masterfile (unzipped) in working directory.
3. Run command with either -m option for coverage metrics only, -v option for
   full variant calling and consensus fasta in offline mode, and new -u option
   to run variant calling and lineage calls using latest version.

``` sh
# Default is offline mode
./run_swift_sarscov2_docker.sh -v sarscov2_v2_masterfile.txt
# Requires online connection
./run_swift_sarscov2_docker.sh -u sarscov2_v2_masterfile.txt
```
