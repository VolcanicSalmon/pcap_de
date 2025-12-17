#tfid.txt fetched from cisbp by 1 organism

while IFS= read -r tf_id
do
  url="https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF=$tf_id"
  echo "Fetching data for $tf_id..."
  curl -s "$url" | grep -A 10 "Position Frequency Matrix" >> tfpf.txt
done <tfid.txt
