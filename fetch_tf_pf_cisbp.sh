#tfid.txt fetched from cisbp by 1 organism

while IFS= read -r tf_id
do
  url="https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF=$tf_id"
  curl -s "$url" | grep -o 'PF[0-9]\+' >> tfpf.txt
done < tfid.txt
