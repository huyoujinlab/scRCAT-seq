# Locate
set -e
FILE_PATH=$(cd "$(dirname "$0")"; pwd)
cd $FILE_PATH/../

# Download data and reference in our pipeline
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=103rREck0vIE6AZvCjPn2Tg_KygDuM78-' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=103rREck0vIE6AZvCjPn2Tg_KygDuM78-" -O data.zip && rm -rf /tmp/cookies.txt
i=`md5sum data.zip | cut -d ' ' -f1`
if [ ! $i = f1523fbd5e0a6a2dabc77b37804b84e9 ]; then
    echo "Error downloading file, please re-run download.sh"
    exit -1
fi
unzip -o data.zip 

# Download genome
mkdir -p ~/index/hg38
cd ~/index/hg38
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
i=`md5sum hg38.fa.gz | cut -d ' ' -f1`
if [ ! $i = 1c9dcaddfa41027f17cd8f7a82c7293b ]; then
    echo "Error downloading genome files, please re-run download.sh"
    exit -1
fi
gzip -d hg38.fa.gz

[[ -d "input" && -d "reference" && -d "script" ]] && echo "right"
[[ -d "input" && -d "reference" && -d "script" ]] || echo "Please check the file, the 'input' folder and 'reference' folder must be in the same directory as the script folde"