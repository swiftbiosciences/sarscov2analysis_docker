# Swift Sarscov2 analysis script

## Pulling image from Docker hub
``` sh
docker pull swiftbiosci/sarscov2analysis:latest
```

## Download Masterfile from Swift Biosciences Bioinformatics Tools and Downloads page

``` sh
# Swift Normalase Amplicon SARS-CoV-2 (SNAP) With additinoal genome coverage Masterfile
curl -O https://swiftbiosci.com/wp-content/uploads/2020/09/sarscov2_v2_masterfile-2.txt.zip
```

## Clone repository swiftbiosciences/sarscov2analysis_docker:

``` sh
git clone https://github.com/swiftbiosciences/sarscov2analysis_docker.git
```

## Wrapper script usage:
1. Copy run_swift_sarscov2_docker.sh to location of your working PATH.
2. Place all fastq files and panel masterfile (unzipped) in working directory.
3. Run command with either -m option for coverage metrics only or -v option
   for full variant calling.
   
``` sh
./run_swift_sarscov2_docker.sh -m sarscov2_v2_masterfile.txt
./run_swift_sarscov2_docker.sh -v sarscov2_v2_masterfile.txt
```
