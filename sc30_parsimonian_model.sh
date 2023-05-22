
ONCHO=/gpfs/gibbs/project/sbsc/ga254/dataproces/ONCHO

for file in $ONCHO/vector_seed*_full/importanceP_selvsVar.txt ; do 
awk '{ if (NR>1) print NR-1 , $1  }'  $file 
done | sort -g -k 1 | uniq -c | awk '{ print $2, $1, $3  }' | sort -g -k 1,1g -k 2,2gr | awk '{ print $2/ $1, $3  }' | sort -k 2,2 > /tmp/test.txt

~/scripts/general/sum.sh /tmp/test.txt  /tmp/importance_all.txt      <<EOF 
n
2
2
EOF

sort -gr -k 1,1 /tmp/importance_all.txt  > $ONCHO/vector/importanceP_all.txt


for COLNAME in  $(awk '{ print $2 }' $ONCHO/vector/importanceP_all.txt) ; do
echo $COLNAME
awk -F, -v COLNAME=$COLNAME ' { if (NR==1)  { for (col=1;col<=NF;col++) { if ($col==COLNAME) {colprint=col}}}  else  {print $1 , $colprint }} ' $ONCHO/vector/x_y_pa_predictors4R.txt    > $ONCHO/vector/predictors4R_selected.txt 
done


