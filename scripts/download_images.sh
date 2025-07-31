# This script is used to download images given a list of urls

imgdir=$1
urllist=$2
idlist=$3

mkdir $imgdir
n=1
while IFS= read -r url && IFS= read -r id <&3;
do
 	wget $url -O $imgdir/$imgdir$id".jpg"
done < $urllist 3< $idlist
echo "Finish downloading"

# deduplicate the urls and remove the headers in the file
# use the command `nohup bash download_images.sh img/ urllist.txt &> download.log &` to run it on the background