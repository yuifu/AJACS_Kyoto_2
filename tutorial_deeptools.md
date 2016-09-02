# AJACS京都２ NGSデータから新たな知識を導出するための高次解析 実習資料


## 実習1: deepTools による aggregation plot の作成

deepTools: http://deeptools.readthedocs.io/en/latest/

```
wget https://www.encodeproject.org/files/ENCFF002CRA/@@download/ENCFF002CRA.bed.gz
gunzip ENCFF002CRA.bed.gz
grep chr21 ENCFF002CRA.bed > ENCFF002CRA.chr21.bed
```

```
wget https://raw.githubusercontent.com/yuifu/AJACS_Kyoto_2/master/tutorial_monocle.md
```


```
bamCoverage -b ENCFF831SAH.chr21.bam -o myfile.bw

computeMatrix reference-point \
       --referencePoint center \
       -b 1000 -a 1000 \
       -R ENCFF002CRA.chr21.bed \
       -S myfile.bw \
       --skipZeros \
       -o matrix_center.gz \
       --outFileSortedRegions regions00.bed

plotProfile -m matrix_center.gz \
              -out ExampleProfile1.png \
              --plotTitle "Test data profile"
```
